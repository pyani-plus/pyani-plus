# The MIT License
#
# Copyright (c) 2024-2025 University of Strathclyde
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""Implements the public user-facing command line interface (CLI).

The commands defined here are intended to be used at the command line, or
perhaps wrapped within a GUI like Galaxy. The CLI should cover running a
new analysis (which can build on existing comparisons if the same DB is
used), and reporting on a finished analysis (exporting tables and plots).
"""

import sys
import tempfile
from contextlib import nullcontext
from math import log
from pathlib import Path
from typing import Annotated

import click
import networkx as nx
import typer
from rich.console import Console
from rich.progress import Progress
from rich.table import Table
from rich.text import Text
from sqlalchemy.orm import Session

from pyani_plus import PROGRESS_BAR_COLUMNS, classify, db_orm, tools
from pyani_plus.methods import anib, anim, fastani, sourmash
from pyani_plus.public_cli_args import (
    OPT_ARG_TYPE_ANIM_MODE,
    OPT_ARG_TYPE_CLASSIFY_MODE,
    OPT_ARG_TYPE_COV_MIN,
    OPT_ARG_TYPE_CREATE_DB,
    OPT_ARG_TYPE_EXECUTOR,
    OPT_ARG_TYPE_FRAGSIZE,
    OPT_ARG_TYPE_KMERSIZE,
    OPT_ARG_TYPE_LABEL,
    OPT_ARG_TYPE_MINMATCH,
    OPT_ARG_TYPE_RUN_ID,
    OPT_ARG_TYPE_RUN_NAME,
    OPT_ARG_TYPE_SOURMASH_SCALED,
    OPT_ARG_TYPE_TEMP,
    OPT_ARG_TYPE_TEMP_WORKFLOW,
    REQ_ARG_TYPE_DATABASE,
    REQ_ARG_TYPE_FASTA_DIR,
    REQ_ARG_TYPE_OUTDIR,
    EnumModeClassify,
)
from pyani_plus.utils import (
    available_cores,
    check_db,
    check_fasta,
    file_md5sum,
    filename_stem,
)
from pyani_plus.workflows import (
    ShowProgress,
    ToolExecutor,
    run_snakemake_with_progress_bar,
)

app = typer.Typer(
    context_settings={"help_option_names": ["-h", "--help"]},
)


def start_and_run_method(  # noqa: PLR0913
    executor: ToolExecutor,
    temp: Path | None,
    workflow_temp: Path | None,
    database: Path,
    name: str | None,
    method: str,
    fasta: Path,
    targets: list[str],  # filenames without paths
    tool: tools.ExternalToolData | None,
    binaries: dict[str, Path],
    *,
    fragsize: int | None = None,
    mode: str | None = None,
    kmersize: int | None = None,
    minmatch: float | None = None,
    extra: str | None = None,
) -> int:
    """Run the snakemake workflow for given method and log run to database."""
    fasta_names = check_fasta(fasta)

    # We should have caught all the obvious failures earlier,
    # including missing inputs or missing external tools.
    # Can now start talking to the DB.
    session = db_orm.connect_to_db(database)

    # Reuse existing config, or log a new one
    config = db_orm.db_configuration(
        session,
        method,
        "" if tool is None else tool.exe_path.stem,
        "" if tool is None else tool.version,
        fragsize,
        mode,
        kmersize,
        minmatch,
        extra,
        create=True,
    )

    n = len(fasta_names)
    filename_to_md5 = {}
    hashes = set()
    with Progress(*PROGRESS_BAR_COLUMNS) as progress:
        for filename in progress.track(fasta_names, description="Indexing FASTAs"):
            md5 = file_md5sum(filename)
            filename_to_md5[filename] = md5
            if md5 in hashes:
                # This avoids hitting IntegrityError UNIQUE constraint failed
                dups = "\n" + "\n".join(
                    sorted({str(k) for k, v in filename_to_md5.items() if v == md5})
                )
                msg = f"ERROR - Multiple genomes with same MD5 checksum {md5}:{dups}"
                sys.exit(msg)
            hashes.add(md5)
            db_orm.db_genome(session, filename, md5, create=True)

    # New run
    run = db_orm.add_run(
        session,
        config,
        cmdline=" ".join(sys.argv),
        fasta_directory=fasta,
        status="Initialising",
        name=f"{len(filename_to_md5)} genomes using {method}" if name is None else name,
        date=None,
        fasta_to_hash=filename_to_md5,
    )
    session.commit()  # Redundant?
    print(f"{method} run setup with {n} genomes in database")

    return run_method(
        executor,
        temp,
        workflow_temp,
        filename_to_md5,
        database,
        session,
        run,
        targets,
        binaries,
    )


def run_method(  # noqa: PLR0913
    executor: ToolExecutor,
    temp: Path | None,
    workflow_temp: Path | None,
    filename_to_md5: dict[Path, str],
    database: Path,
    session: Session,
    run: db_orm.Run,
    targets: list[str],  # filenames without paths
    binaries: dict[str, Path],
) -> int:
    """Run the snakemake workflow for given method and log run to database."""
    run_id = run.run_id
    configuration = run.configuration
    method = configuration.method
    workflow_name = (
        f"snakemake_{method.lower()}.smk"
        if method in ("sourmash",)
        else "compute_column.smk"
    )
    params: dict[str, object] = {
        # Paths etc - see also outdir below
        "indir": Path(run.fasta_directory).resolve(),  # must be absolute
        "db": Path(database).resolve(),  # must be absolute
        "run_id": run_id,
        "cores": available_cores(),  # should make configurable
        # Method settings:
        "fragsize": configuration.fragsize,
        "mode": configuration.mode,
        "kmersize": configuration.kmersize,
        "minmatch": configuration.minmatch,
        "extra": configuration.extra,
        "md5_to_filename": {k: str(v) for (v, k) in filename_to_md5.items()},
    }
    params.update({k: str(v) for k, v in binaries.items()})
    del configuration

    done = run.comparisons().count()
    n = len(filename_to_md5)
    if done == n**2:
        print(f"Database already has all {n}²={n**2} {method} comparisons")
    else:
        print(
            f"Database already has {done} of {n}²={n**2} {method} comparisons, {n**2 - done} needed"
        )
        run.status = "Running"
        session.commit()
        session.close()  # Reduce chance of DB locking
        del run

        # Run snakemake wrapper
        # With a cluster-based running like SLURM, the location of the working and
        # output directories must be viewable from all the worker nodes too.
        # i.e. Can't use a temp directory on the head node.
        # We might want to make this explicitly configurable, e.g. to /mnt/scratch/
        with (
            nullcontext(workflow_temp)
            if workflow_temp
            else tempfile.TemporaryDirectory(
                prefix="pyani-plus_", dir=None if executor.value == "local" else "."
            ) as tmp
        ):
            work_path = Path(tmp) / "working"
            out_path = Path(tmp) / "output"
            params["outdir"] = out_path.resolve()
            target_paths = [out_path / _ for _ in targets]
            run_snakemake_with_progress_bar(
                executor,
                workflow_name,
                target_paths,
                params,
                work_path,
                display=ShowProgress.spin if method == "sourmash" else ShowProgress.bar,
                database=Path(database),
                run_id=run_id,
                temp=temp,
            )

        # Reconnect to the DB
        session = db_orm.connect_to_db(database)
        run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()
        done = run.comparisons().count()

    if done != n**2:
        # There is no obvious way to test this hypothetical failure:
        msg = f"ERROR: Only have {done} of {n}²={n**2} {method} comparisons needed"  # pragma: no cover
        sys.exit(msg)  # pragma: no cover

    run.cache_comparisons()  # will this needs a progress bar too with large n?
    run.status = "Done"
    session.commit()
    print(f"Completed {method} run-id {run_id} with {n} genomes in database {database}")
    session.close()
    return 0


@app.command("anim", rich_help_panel="ANI methods")
def cli_anim(  # noqa: PLR0913
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    *,
    # These are for the run table:
    name: OPT_ARG_TYPE_RUN_NAME = None,
    # Does not use fragsize, kmersize, or minmatch
    # The mode here is not optional - must pick one!
    mode: OPT_ARG_TYPE_ANIM_MODE = anim.MODE,
    create_db: OPT_ARG_TYPE_CREATE_DB = False,
    executor: OPT_ARG_TYPE_EXECUTOR = ToolExecutor.local,
    temp: OPT_ARG_TYPE_TEMP = None,
    wtemp: OPT_ARG_TYPE_TEMP_WORKFLOW = None,
) -> int:
    """Execute ANIm calculations, logged to a pyANI-plus SQLite3 database."""
    check_db(database, create_db)

    tool = tools.get_nucmer()
    binaries = {
        "nucmer": tool.exe_path,
        "delta_filter": tools.get_delta_filter().exe_path,
    }
    fasta_list = check_fasta(fasta)

    return start_and_run_method(
        executor,
        temp,
        wtemp,
        database,
        name,
        "ANIm",
        fasta,
        [f"all_vs_{Path(_).stem}.anim" for _ in fasta_list],
        tool,
        binaries,
        mode=mode.value,  # turn the enum into a string
    )


@app.command("dnadiff", rich_help_panel="ANI methods")
def cli_dnadiff(  # noqa: PLR0913
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    *,
    # These are for the run table:
    name: OPT_ARG_TYPE_RUN_NAME = None,
    # Does not use fragsize, mode, kmersize, or minmatch
    create_db: OPT_ARG_TYPE_CREATE_DB = False,
    executor: OPT_ARG_TYPE_EXECUTOR = ToolExecutor.local,
    temp: OPT_ARG_TYPE_TEMP = None,
    wtemp: OPT_ARG_TYPE_TEMP_WORKFLOW = None,
) -> int:
    """Execute mumer-based dnadiff calculations, logged to a pyANI-plus SQLite3 database."""
    check_db(database, create_db)

    # We don't actually call the tool dnadiff (which has its own version),
    # rather we call nucmer, delta-filter, show-diff and show-coords from MUMmer
    tool = tools.get_nucmer()
    binaries = {
        "nucmer": tool.exe_path,
        "delta_filter": tools.get_delta_filter().exe_path,
        "show_diff": tools.get_show_diff().exe_path,
        "show_coords": tools.get_show_coords().exe_path,
    }
    fasta_list = check_fasta(fasta)

    return start_and_run_method(
        executor,
        temp,
        wtemp,
        database,
        name,
        "dnadiff",
        fasta,
        [f"all_vs_{Path(_).stem}.dnadiff" for _ in fasta_list],
        tool,
        binaries,
    )


@app.command("anib", rich_help_panel="ANI methods")
def cli_anib(  # noqa: PLR0913
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the run table:
    *,
    name: OPT_ARG_TYPE_RUN_NAME = None,
    # These are all for the configuration table:
    fragsize: OPT_ARG_TYPE_FRAGSIZE = anib.FRAGSIZE,
    # Does not use mode, kmersize, or minmatch
    create_db: OPT_ARG_TYPE_CREATE_DB = False,
    executor: OPT_ARG_TYPE_EXECUTOR = ToolExecutor.local,
    temp: OPT_ARG_TYPE_TEMP = None,
    wtemp: OPT_ARG_TYPE_TEMP_WORKFLOW = None,
) -> int:
    """Execute ANIb calculations, logged to a pyANI-plus SQLite3 database."""
    check_db(database, create_db)

    tool = tools.get_blastn()
    alt = tools.get_makeblastdb()
    if tool.version != alt.version:
        msg = f"ERROR: blastn {tool.version} vs makeblastdb {alt.version}"
        sys.exit(msg)
    binaries = {
        "blastn": tool.exe_path,
        "makeblastdb": alt.exe_path,
    }
    fasta_list = check_fasta(fasta)

    return start_and_run_method(
        executor,
        temp,
        wtemp,
        database,
        name,
        "ANIb",
        fasta,
        [f"all_vs_{Path(_).stem}.anib" for _ in fasta_list],
        tool,
        binaries,
        fragsize=fragsize,
    )


@app.command("fastani", rich_help_panel="ANI methods")
def cli_fastani(  # noqa: PLR0913
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    *,
    # These are for the run table:
    name: OPT_ARG_TYPE_RUN_NAME = None,
    # These are all for the configuration table:
    fragsize: OPT_ARG_TYPE_FRAGSIZE = fastani.FRAG_LEN,
    # Does not use mode
    # Don't use OPT_ARG_TYPE_KMERSIZE as want to include max=16
    kmersize: Annotated[
        int,
        typer.Option(
            help="Comparison method k-mer size",
            rich_help_panel="Method parameters",
            min=1,
            max=16,
        ),
    ] = fastani.KMER_SIZE,
    minmatch: OPT_ARG_TYPE_MINMATCH = fastani.MIN_FRACTION,
    create_db: OPT_ARG_TYPE_CREATE_DB = False,
    executor: OPT_ARG_TYPE_EXECUTOR = ToolExecutor.local,
    temp: OPT_ARG_TYPE_TEMP = None,
    wtemp: OPT_ARG_TYPE_TEMP_WORKFLOW = None,
) -> int:
    """Execute fastANI calculations, logged to a pyANI-plus SQLite3 database."""
    check_db(database, create_db)

    tool = tools.get_fastani()
    binaries = {
        "fastani": tool.exe_path,
    }
    fasta_list = check_fasta(fasta)

    return start_and_run_method(
        executor,
        temp,
        wtemp,
        database,
        name,
        "fastANI",
        fasta,
        [f"all_vs_{Path(_).stem}.fastani" for _ in fasta_list],
        tool,
        binaries,
        fragsize=fragsize,
        kmersize=kmersize,
        minmatch=minmatch,
    )


@app.command("sourmash", rich_help_panel="ANI methods")
def cli_sourmash(  # noqa: PLR0913
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    *,
    # These are for the run table:
    name: OPT_ARG_TYPE_RUN_NAME = None,
    create_db: OPT_ARG_TYPE_CREATE_DB = False,
    executor: OPT_ARG_TYPE_EXECUTOR = ToolExecutor.local,
    temp: OPT_ARG_TYPE_TEMP = None,
    wtemp: OPT_ARG_TYPE_TEMP_WORKFLOW = None,
    # These are for the configuration table:
    scaled: OPT_ARG_TYPE_SOURMASH_SCALED = sourmash.SCALED,  # 1000
    kmersize: OPT_ARG_TYPE_KMERSIZE = sourmash.KMER_SIZE,
) -> int:
    """Execute sourmash-plugin-branchwater ANI calculations, logged to a pyANI-plus SQLite3 database."""
    check_db(database, create_db)

    tool = tools.get_sourmash()
    binaries = {
        "sourmash": tool.exe_path,
    }
    extra = f"scaled={scaled}"
    return start_and_run_method(
        executor,
        temp,
        wtemp,
        database,
        name,
        "sourmash",
        fasta,
        ["manysearch.csv"],
        tool,
        binaries,
        kmersize=kmersize,
        extra=extra,
    )


@app.command(rich_help_panel="ANI methods")
def external_alignment(  # noqa: PLR0913
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    *,
    # These are for the run table:
    name: OPT_ARG_TYPE_RUN_NAME = None,
    create_db: OPT_ARG_TYPE_CREATE_DB = False,
    executor: OPT_ARG_TYPE_EXECUTOR = ToolExecutor.local,
    temp: OPT_ARG_TYPE_TEMP = None,
    wtemp: OPT_ARG_TYPE_TEMP_WORKFLOW = None,
    # These are for the configuration table:
    alignment: Annotated[
        Path,
        typer.Option(
            help="FASTA format MSA of the same genomes (one sequence per genome)",
            rich_help_panel="Method parameters",
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ],
    label: Annotated[
        str,
        typer.Option(
            click_type=click.Choice(["md5", "filename", "stem"]),
            rich_help_panel="Method parameters",
            help="How are the sequences in the MSA labelled vs the FASTA genomes?",
        ),
    ] = "stem",
) -> int:
    """Compute pairwise ANI from given multiple-sequence-alignment (MSA) file."""
    check_db(database, create_db)

    aln_checksum = file_md5sum(alignment)
    # Doing this order to put the filename LAST, in case of separators in the filename
    extra = f"md5={aln_checksum};label={label};alignment={alignment.name}"
    fasta_list = check_fasta(fasta)

    return start_and_run_method(
        executor,
        temp,  # not needed?
        wtemp,  # not needed?
        database,
        f"Import of {alignment.name}" if name is None else name,
        "external-alignment",
        fasta,
        [f"all_vs_{Path(_).stem}.external-alignment" for _ in fasta_list],
        None,  # no tool
        {},  # no binaries
        extra=extra,
    )


@app.command()
def resume(  # noqa: C901, PLR0912, PLR0915
    database: REQ_ARG_TYPE_DATABASE,
    *,
    run_id: OPT_ARG_TYPE_RUN_ID = None,
    executor: OPT_ARG_TYPE_EXECUTOR = ToolExecutor.local,
    temp: OPT_ARG_TYPE_TEMP = None,
    wtemp: OPT_ARG_TYPE_TEMP_WORKFLOW = None,
) -> int:
    """Resume any (partial) run already logged in the database.

    If the run was already complete, this should have no effect.

    Any missing pairwise comparisons will be computed, and the the old
    run will be marked as complete.

    If the version of the underlying tool has changed, this will abort
    as the original run cannot be completed.
    """
    if database == ":memory:" or not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    session = db_orm.connect_to_db(database)
    run = db_orm.load_run(session, run_id)
    if run_id is None:
        run_id = run.run_id  # relevant if was None
        print(f"INFO: Resuming run-id {run_id}")
    config = run.configuration
    print(
        f"INFO: This is a {config.method} run on {run.genomes.count()} genomes, "
        f"using {config.program} version {config.version}"
    )
    if not run.genomes.count():
        msg = f"ERROR: No genomes recorded for run-id {run_id}, cannot resume."
        sys.exit(msg)

    # The params dict has two kinds of entries,
    # - tool paths, which ought to be handled more neatly
    # - config entries, which ought to be named consistently and done centrally
    tool: tools.ExternalToolData | None = None  # make type explicit for mypy
    match config.method:
        case "fastANI":
            tool = tools.get_fastani()
            binaries = {
                "fastani": tool.exe_path,
            }
            # Ideally this would be selective, as this will recompute every column
            targets = [
                f"all_vs_{Path(_.fasta_filename).stem}.fastani"
                for _ in run.fasta_hashes
            ]
        case "ANIm":
            tool = tools.get_nucmer()
            binaries = {
                "nucmer": tool.exe_path,
                "delta_filter": tools.get_delta_filter().exe_path,
            }
            targets = [
                f"all_vs_{Path(_.fasta_filename).stem}.anim" for _ in run.fasta_hashes
            ]
        case "dnadiff":
            tool = tools.get_nucmer()
            binaries = {
                "nucmer": tool.exe_path,
                "delta_filter": tools.get_delta_filter().exe_path,
                "show_diff": tools.get_show_diff().exe_path,
                "show_coords": tools.get_show_coords().exe_path,
            }
            targets = [
                f"all_vs_{Path(_.fasta_filename).stem}.dnadiff"
                for _ in run.fasta_hashes
            ]
        case "ANIb":
            tool = tools.get_blastn()
            binaries = {
                "blastn": tool.exe_path,
                "makeblastdb": tools.get_makeblastdb().exe_path,
            }
            targets = [
                f"all_vs_{Path(_.fasta_filename).stem}.anib" for _ in run.fasta_hashes
            ]
        case "sourmash":
            tool = tools.get_sourmash()
            binaries = {
                "sourmash": tool.exe_path,
            }
            targets = ["manysearch.csv"]
        case "external-alignment":
            tool = None
            binaries = {}
            targets = [
                f"all_vs_{Path(_.fasta_filename).stem}.external-alignment"
                for _ in run.fasta_hashes
            ]
        case _:
            msg = f"ERROR: Unknown method {config.method} for run-id {run_id} in {database}"
            sys.exit(msg)
    if not tool:
        if config.program != "" or config.version != "":
            msg = (
                "ERROR: We expect no tool information, but"
                f" run-id {run_id} used {config.program} version {config.version} instead."
            )
            sys.exit(msg)
    elif tool.exe_path.stem != config.program or tool.version != config.version:
        msg = (
            f"ERROR: We have {tool.exe_path.stem} version {tool.version}, but"
            f" run-id {run_id} used {config.program} version {config.version} instead."
        )
        sys.exit(msg)

    del tool

    # Now we need to check the fasta files in the directory
    # against those included in the run...
    fasta = Path(run.fasta_directory)
    if not fasta.is_dir():
        msg = f"ERROR: run-id {run_id} used input folder {fasta}, but that is not a directory (now)."
        sys.exit(msg)

    # Recombine the fasta directory name from the runs table with the plain filename from
    # the run-genome linking table
    filename_to_md5 = {
        fasta / link.fasta_filename: link.genome_hash for link in run.fasta_hashes
    }
    for filename, md5 in filename_to_md5.items():
        if not filename.is_file():
            msg = (
                f"ERROR: run-id {run_id} used {filename} with MD5 {md5}"
                f" but this FASTA file no longer exists"
            )
            sys.exit(msg)

    # Resuming
    run.status = "Resuming"
    session.commit()

    return run_method(
        executor,
        temp,
        wtemp,
        filename_to_md5,
        database,
        session,
        run,
        targets,
        binaries,
    )


@app.command()
def list_runs(
    database: REQ_ARG_TYPE_DATABASE,
) -> int:
    """List the runs defined in a given pyANI-plus SQLite3 database."""
    if database == ":memory:" or not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    session = db_orm.connect_to_db(database)
    runs = session.query(db_orm.Run)

    table = Table(
        title=f"{runs.count()} analysis runs in {database}",
        row_styles=["dim", ""],  # alternating zebra stripes
    )
    table.add_column("ID", justify="right", no_wrap=True)
    table.add_column("Date")
    table.add_column("Method")
    table.add_column(Text("Done", justify="left"), justify="right", no_wrap=True)
    table.add_column(Text("Null", justify="left"), justify="right", no_wrap=True)
    table.add_column(Text("Miss", justify="left"), justify="right", no_wrap=True)
    table.add_column(Text("Total", justify="left"), justify="right", no_wrap=True)
    table.add_column("Status")
    table.add_column("Name")
    # Would be nice to report {conf.program} {conf.version} too,
    # perhaps conditional on the terminal width?
    for run in runs:
        conf = run.configuration
        n = run.genomes.count()
        total = n**2
        done = run.comparisons().count()
        # Using is None does not work as expected, must use == None
        nulls = run.comparisons().where(db_orm.Comparison.identity == None).count()  # noqa: E711
        table.add_row(
            str(run.run_id),
            str(run.date.date()),
            conf.method,
            Text(
                f"{done - nulls}",
                style="green" if done == total and not nulls else "yellow",
            ),
            Text(f"{nulls}", style="red" if nulls else "green"),
            Text(f"{total - done}", style="yellow" if done < total else "green"),
            f"{n**2}={n}²",
            run.status,
            run.name,
        )
    session.close()
    console = Console()
    console.print(table)
    return 0


@app.command()
def delete_run(
    database: REQ_ARG_TYPE_DATABASE,
    run_id: OPT_ARG_TYPE_RUN_ID = None,
    *,
    force: Annotated[
        # Listing name(s) explicitly to avoid automatic matching --no-create-db
        bool, typer.Option("-f", "--force", help="Delete without confirmation")
    ] = False,
) -> int:
    """Delete any single run from the given pyANI-plus SQLite3 database.

    This will prompt the user for confirmation if the run has comparisons, or
    if the run status is "Running", but that can be overridden.

    Currently this will *not* delete any linked comparisons, even if they are
    not currently linked to another run. They will be reused should you start
    a new run using an overlapping set of input FASTA files.
    """
    if database == ":memory:" or not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    confirm = False

    session = db_orm.connect_to_db(database)
    run = db_orm.load_run(session, run_id, check_complete=False)
    if run_id is None:
        run_id = run.run_id
        print("INFO: Deleting most recent run")
        confirm = True

    # Could use rish colours, match how list-runs colours the counts?
    done = run.comparisons().count()
    n = run.genomes.count()
    if n and done == n**2:
        print(
            f"INFO: Run {run_id} contains all {n**2}={n}²"
            f" {run.configuration.method} comparisons, status: {run.status}"
        )
        confirm = True
    else:
        print(
            f"INFO: Run {run_id} contains {done}/{n**2}={n}²"
            f" {run.configuration.method} comparisons, status: {run.status}"
        )
        if done:
            confirm = True
    print(f"INFO: Run name: {run.name}")

    if run.status == "Running":  # should be a constant or enum?
        # Should we also look at the date of the run? If old probably it failed.
        print("WARNING: Deleting a run still being computed will cause it to fail!")
        confirm = True

    if confirm and not force:
        click.confirm("Do you want to continue?", abort=True)  # pragma: no cover

    # Plan to offer an extended mode, perhaps --all, which will also
    # delete orphaned entries in the configuration,  genomes, and
    # comparisons tables. By default will just drop the runs entry
    # and direct links in runs_genomes thanks to the relationship.

    # Explicitly delete the no longer wanted entries in runs_genomes
    # (ought to be able to trigger this automatically in a cascade
    # when the run itself is deleted, but that didn't work for me):
    session.query(db_orm.RunGenomeAssociation).where(
        db_orm.RunGenomeAssociation.run_id == run_id
    ).delete()
    session.delete(run)
    session.commit()
    session.close()
    return 0


@app.command()
def export_run(  # noqa: C901
    database: REQ_ARG_TYPE_DATABASE,
    outdir: REQ_ARG_TYPE_OUTDIR,
    run_id: OPT_ARG_TYPE_RUN_ID = None,
    label: OPT_ARG_TYPE_LABEL = "stem",
) -> int:
    """Export any single run from the given pyANI-plus SQLite3 database.

    The output directory must already exist. Any pre-existing files will be
    overwritten.

    The matrix output files are named <method>_<property>.tsv while the
    long form is named <method>_run_<run-id>.tsv and will include the
    query and subject genomes and all the comparison properties as columns.

    Incomplete runs will return an error. There will be no output for empty
    run. For partial runs the long form table will be exported, but not the
    matrices.
    """
    if not outdir.is_dir():
        msg = f"ERROR: Output directory {outdir} does not exist"
        sys.exit(msg)

    if database == ":memory:" or not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    session = db_orm.connect_to_db(database)
    run = db_orm.load_run(session, run_id, check_empty=True)
    if run_id is None:
        run_id = run.run_id
        print(f"INFO: Exporting run-id {run_id}")

    # Question: Should we export a plain text of JSON summary of the configuration etc?
    # Question: Should we include the run-id in the matrix filenames?
    # Question: Should we match the property in filenames to old pyANI (e.g. coverage)?
    method = run.configuration.method

    def float_or_na(value: float | None) -> str:
        """Format a float where represent None as NA (null in our DB)."""
        # Should this allow configurable float formatting?
        return "NA" if value is None else str(value)

    if label == "md5":
        mapping = lambda x: x  # noqa: E731
    elif label == "filename":
        mapping = {_.genome_hash: _.fasta_filename for _ in run.fasta_hashes}.get
    else:
        mapping = {
            _.genome_hash: filename_stem(_.fasta_filename) for _ in run.fasta_hashes
        }.get

    long_filename = f"{method}_run_{run_id}.tsv"
    with (outdir / long_filename).open("w") as handle:
        # Should the column names match our internal naming?
        handle.write(
            "#Query\tSubject\tIdentity\tQuery-Cov\tSubject-Cov\tHadamard\ttANI\tAlign-Len\tSim-Errors\n"
        )
        for _ in run.comparisons():
            # Below for the matrix output we use the cached DataFrame for Hadamard.
            # Here compute it on the fly (we might be exporting partial results).
            hadamard = (
                None
                if _.identity is None or _.cov_query is None
                else _.identity * _.cov_query
            )
            tani = None if hadamard is None else -log(hadamard)
            handle.write(
                f"{mapping(_.query_hash)}\t{mapping(_.subject_hash)}"
                f"\t{float_or_na(_.identity)}"
                f"\t{float_or_na(_.cov_query)}"
                f"\t{float_or_na(_.cov_subject)}"
                f"\t{float_or_na(hadamard)}"
                f"\t{float_or_na(tani)}"
                f"\t{float_or_na(_.aln_length)}"
                f"\t{float_or_na(_.sim_errors)}\n"
            )
    print(f"Wrote long-form to {outdir}/{long_filename}")

    # Reload the run checking it is complete (quick) (might abort here!),
    # and caching the matrices if needed (slower):
    run = db_orm.load_run(session, run_id, check_complete=True)

    for matrix, filename in (
        (run.identities, f"{method}_identity.tsv"),
        (run.aln_length, f"{method}_aln_lengths.tsv"),
        (run.sim_errors, f"{method}_sim_errors.tsv"),
        (run.cov_query, f"{method}_query_cov.tsv"),
        (run.hadamard, f"{method}_hadamard.tsv"),
        (run.tani, f"{method}_tANI.tsv"),
    ):
        if matrix is None:
            # This is mainly for mypy to assert the matrix is not None
            msg = f"ERROR: Could not load run {method} matrix"  # pragma: no cover
            sys.exit(msg)  # pragma: no cover

        try:
            matrix = run.relabelled_matrix(matrix, label)  # noqa: PLW2901
        except ValueError as err:
            msg = f"ERROR: {err}"
            sys.exit(msg)

        matrix.to_csv(outdir / filename, sep="\t")

    print(f"Wrote matrices to {outdir}/{method}_*.tsv")

    session.close()
    return 0


@app.command()
def plot_run(
    database: REQ_ARG_TYPE_DATABASE,
    outdir: REQ_ARG_TYPE_OUTDIR,
    run_id: OPT_ARG_TYPE_RUN_ID = None,
    label: OPT_ARG_TYPE_LABEL = "stem",
) -> int:
    """Plot heatmaps and distributions for any single run.

    The output directory must already exist. The heatmap files will be named
    <method>_<property>.<extension> and any pre-existing files will be overwritten.
    """
    if not outdir.is_dir():
        msg = f"ERROR: Output directory {outdir} does not exist"
        sys.exit(msg)

    if database == ":memory:" or not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    session = db_orm.connect_to_db(database)
    run = db_orm.load_run(session, run_id, check_complete=True)
    if run_id is None:
        run_id = run.run_id
        print(f"INFO: Plotting {run.configuration.method} run-id {run_id}")

    from pyani_plus import plot_run  # lazy import

    count = plot_run.plot_single_run(run, outdir, label)
    print(f"Wrote {count} images to {outdir}/{run.configuration.method}_*.*")
    session.close()
    return 0 if count else 1


@app.command()
def plot_run_comp(
    database: REQ_ARG_TYPE_DATABASE,
    outdir: REQ_ARG_TYPE_OUTDIR,
    run_ids: Annotated[
        str,
        typer.Option(help="Which runs (comma separated list, reference first)?"),
    ],
    columns: Annotated[
        int,
        typer.Option(
            help="How many columns to use when tiling plots of multiple runs."
            " Default 0 means automatically tries for square tiling.",
            min=0,
        ),
    ] = 0,
) -> int:
    """Plot comparisons between multiple runs.

    The output directory must already exist. The scatter plots will be named
    <method>_<property>_<run-id>_vs_*.<extension> and any
    pre-existing files will be overwritten.
    """
    if not outdir.is_dir():
        msg = f"ERROR: Output directory {outdir} does not exist"
        sys.exit(msg)

    if database == ":memory:" or not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    try:
        runs = [int(_) for _ in run_ids.split(",")]
    except ValueError:
        msg = f"ERROR: Expected comma separated list of runs, not: {run_ids}"
        sys.exit(msg)

    run_id = runs[0]  # the reference
    other_runs = runs[1:]

    if not other_runs:
        msg = "ERROR: Need at least two runs for a comparison"
        sys.exit(msg)

    session = db_orm.connect_to_db(database)
    ref_run = db_orm.load_run(session, run_id, check_complete=False)

    if not ref_run.comparisons().count():
        msg = f"ERROR: Run {run_id} has no comparisons"
        sys.exit(msg)

    from pyani_plus import plot_run  # lazy import

    done = plot_run.plot_run_comparison(session, ref_run, other_runs, outdir, columns)
    print(
        f"Wrote {done} images to {outdir}/{ref_run.configuration.method}_identity_{run_id}_vs_*.*"
    )
    session.close()
    return 0


@app.command("classify", rich_help_panel="Commands")
def cli_classify(  # noqa: C901, PLR0912, PLR0913, PLR0915
    database: REQ_ARG_TYPE_DATABASE,
    outdir: REQ_ARG_TYPE_OUTDIR,
    coverage_edges: Annotated[
        str,
        typer.Option(
            help="How to resolve asymmetrical ANI coverage results for edges in the graph (min, max or mean).",
            rich_help_panel="Method parameters",
        ),
    ] = "min",
    score_edges: Annotated[
        str,
        typer.Option(
            help="How to resolve asymmetrical ANI identity/tANI results for edges in the graph (min, max or mean).",
            rich_help_panel="Method parameters",
        ),
    ] = "mean",
    run_id: OPT_ARG_TYPE_RUN_ID = None,
    label: OPT_ARG_TYPE_LABEL = "stem",
    cov_min: OPT_ARG_TYPE_COV_MIN = classify.MIN_COVERAGE,
    mode: OPT_ARG_TYPE_CLASSIFY_MODE = classify.MODE,
) -> int:
    """Classify genomes into clusters based on ANI results."""
    if not outdir.is_dir():
        msg = f"ERROR: Output directory {outdir} does not exist"
        sys.exit(msg)

    if database == ":memory:" or not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    session = db_orm.connect_to_db(database)
    run = db_orm.load_run(session, run_id, check_complete=True)
    if run_id is None:
        run_id = run.run_id
        print(f"INFO: Exporting run-id {run_id}")

    method = run.configuration.method

    matrix = None
    if mode == "identity":
        matrix = run.identities
    elif mode == "tANI":
        matrix = run.tani

    if matrix is None:
        msg = f"ERROR: Could not load run {method} matrix"  # pragma: no cover
        sys.exit(msg)  # pragma: no cover

    done = run.comparisons().count()
    run_genomes = run.genomes.count()

    single_genome_run = False
    if done == 1 and run_genomes == 1:
        msg = f"WARNING: Run {run_id} has {done} comparison across {run_genomes} genome. Reporting single clique...\n"
        single_genome_run = True
        sys.stderr.write(msg)  # pragma: no cover
    else:
        print(
            f"Run {run_id} has {done} comparisons across {run_genomes} genomes. Running classify..."
        )

    cov = run.cov_query
    if cov is None:
        msg = f"ERROR: Could not load run {method} matrix"  # pragma: no cover
        sys.exit(msg)  # pragma: no cover

    try:
        score_matrix = run.relabelled_matrix(matrix, label)
        cov = run.relabelled_matrix(cov, label)
    except ValueError as err:
        msg = f"ERROR: {err}"
        sys.exit(msg)

    # Map the string inputs to callable functions
    coverage_agg_func = classify.AGG_FUNCS[coverage_edges]
    identity_agg_func = classify.AGG_FUNCS[score_edges]

    # Construct the graph with the correct functions
    complete_graph = classify.construct_graph(
        cov, score_matrix, coverage_agg_func, identity_agg_func, cov_min
    )
    # Finding cliques
    if len(list(nx.connected_components(complete_graph))) != 1:
        initial_cliques = classify.find_initial_cliques(complete_graph)
    else:
        initial_cliques = []
    recursive_cliques = classify.find_cliques_recursively(complete_graph)

    # Get a list of unique cliques to avoid duplicates. Prioritise initial_cliques
    unique_cliques = classify.get_unique_cliques(initial_cliques, recursive_cliques)

    # Determine column name based on mode
    suffix = "identity" if mode == EnumModeClassify.identity else "tANI"
    column_map = {
        "min_score": f"min_{suffix}",
        "max_score": f"max_{suffix}",
    }

    # Writing the results to .tsv
    clique_data, clique_df = classify.compute_classify_output(
        unique_cliques, method, outdir, column_map
    )
    print(f"Wrote classify output to {outdir}")

    # Only plot classify if more than one genome in comparisons
    # NOTE TO SELF: This will still break if all genomes are siglentons!
    if not single_genome_run:
        genome_groups = classify.get_genome_cligue_ids(clique_df, suffix)
        genome_positions = classify.get_genome_order(genome_groups)
        classify.plot_classify(genome_positions, clique_df, outdir, method, suffix)
    session.close()
    return 0


if __name__ == "__main__":
    sys.exit(app())  # pragma: no cover
