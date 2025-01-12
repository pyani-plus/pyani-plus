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
from pathlib import Path
from typing import Annotated

import click
import typer
from rich.console import Console
from rich.progress import Progress
from rich.table import Table
from rich.text import Text
from sqlalchemy.orm import Session

from pyani_plus import PROGRESS_BAR_COLUMNS, db_orm, tools
from pyani_plus.methods import method_anib, method_anim, method_fastani, method_sourmash
from pyani_plus.public_cli_args import (
    OPT_ARG_TYPE_ANIM_MODE,
    OPT_ARG_TYPE_CREATE_DB,
    OPT_ARG_TYPE_EXECUTOR,
    OPT_ARG_TYPE_FRAGSIZE,
    OPT_ARG_TYPE_KMERSIZE,
    OPT_ARG_TYPE_MINMATCH,
    OPT_ARG_TYPE_RUN_ID,
    OPT_ARG_TYPE_RUN_NAME,
    OPT_ARG_TYPE_SOURMASH_SCALED,
    OPT_ARG_TYPE_TEMP,
    OPT_ARG_TYPE_TEMP_WORKFLOW,
    REQ_ARG_TYPE_DATABASE,
    REQ_ARG_TYPE_FASTA_DIR,
    REQ_ARG_TYPE_OUTDIR,
)
from pyani_plus.utils import (
    available_cores,
    check_db,
    check_fasta,
    file_md5sum,
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
    tool: tools.ExternalToolData,
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
        tool.exe_path.stem,
        tool.version,
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
    workflow_name = f"snakemake_{method.lower()}.smk"
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
                display=ShowProgress.spin
                if method in {"sourmash", "branchwater"}
                else ShowProgress.bar,
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


@app.command(rich_help_panel="ANI methods")
def anim(  # noqa: PLR0913
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    *,
    # These are for the run table:
    name: OPT_ARG_TYPE_RUN_NAME = None,
    # Does not use fragsize, kmersize, or minmatch
    # The mode here is not optional - must pick one!
    mode: OPT_ARG_TYPE_ANIM_MODE = method_anim.MODE,
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


@app.command(rich_help_panel="ANI methods")
def dnadiff(  # noqa: PLR0913
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


@app.command(rich_help_panel="ANI methods")
def anib(  # noqa: PLR0913
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the run table:
    *,
    name: OPT_ARG_TYPE_RUN_NAME = None,
    # These are all for the configuration table:
    fragsize: OPT_ARG_TYPE_FRAGSIZE = method_anib.FRAGSIZE,
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


@app.command(rich_help_panel="ANI methods")
def fastani(  # noqa: PLR0913
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    *,
    # These are for the run table:
    name: OPT_ARG_TYPE_RUN_NAME = None,
    # These are all for the configuration table:
    fragsize: OPT_ARG_TYPE_FRAGSIZE = method_fastani.FRAG_LEN,
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
    ] = method_fastani.KMER_SIZE,
    minmatch: OPT_ARG_TYPE_MINMATCH = method_fastani.MIN_FRACTION,
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


@app.command(rich_help_panel="ANI methods")
def sourmash(  # noqa: PLR0913
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    *,
    # These are for the run table:
    name: OPT_ARG_TYPE_RUN_NAME = None,
    create_db: OPT_ARG_TYPE_CREATE_DB = False,
    executor: OPT_ARG_TYPE_EXECUTOR = ToolExecutor.local,
    temp: OPT_ARG_TYPE_TEMP = None,
    wtemp: OPT_ARG_TYPE_TEMP_WORKFLOW = None,
    # For the config table:
    scaled: OPT_ARG_TYPE_SOURMASH_SCALED = method_sourmash.SCALED,  # 1000
    kmersize: OPT_ARG_TYPE_KMERSIZE = method_sourmash.KMER_SIZE,
) -> int:
    """Execute sourmash calculations, logged to a pyANI-plus SQLite3 database."""
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
        # Can we do e.g. sourmash_k=31_scaled=300.csv
        # using [f"sourmash_k={kmersize}_{extra}.csv"] ?
        ["sourmash.csv"],
        tool,
        binaries,
        kmersize=kmersize,
        extra=extra,
    )


@app.command(rich_help_panel="ANI methods")
def branchwater(  # noqa: PLR0913
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
    scaled: OPT_ARG_TYPE_SOURMASH_SCALED = method_sourmash.SCALED,  # 1000
    kmersize: OPT_ARG_TYPE_KMERSIZE = method_sourmash.KMER_SIZE,
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
        "branchwater",
        fasta,
        ["manysearch.csv"],
        tool,
        binaries,
        kmersize=kmersize,
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
            targets = ["sourmash.csv"]
        case "branchwater":
            tool = tools.get_sourmash()
            binaries = {
                "sourmash": tool.exe_path,
            }
            targets = ["manysearch.csv"]
        case _:
            msg = f"ERROR: Unknown method {config.method} for run-id {run_id} in {database}"
            sys.exit(msg)
    if tool.exe_path.stem != config.program or tool.version != config.version:
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
def export_run(
    database: REQ_ARG_TYPE_DATABASE,
    outdir: REQ_ARG_TYPE_OUTDIR,
    run_id: OPT_ARG_TYPE_RUN_ID = None,
    # Would like to replace this with Literal["md5", "filename", "stem"] once typer updated
    label: Annotated[
        str,
        typer.Option(
            click_type=click.Choice(["md5", "filename", "stem"]),
            help="How to label the genomes",
        ),
    ] = "stem",
) -> int:
    """Export any single run from the given pyANI-plus SQLite3 database.

    The output directory must already exist. The output files are named
    <method>_<property>.tsv and any pre-existing files will be overwritten.
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
        print(f"INFO: Exporting run-id {run_id}")

    # Question: Should we export a plain text of JSON summary of the configuration etc?
    # Question: Should we include the run-id in the filenames name?
    # Question: Should we match the property in filenames to old pyANI (e.g. coverage)?
    method = run.configuration.method

    for matrix, filename in (
        (run.identities, f"{method}_identity.tsv"),
        (run.aln_length, f"{method}_aln_lengths.tsv"),
        (run.sim_errors, f"{method}_sim_errors.tsv"),
        (run.cov_query, f"{method}_query_cov.tsv"),
        (run.hadamard, f"{method}_hadamard.tsv"),
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
def plot_run(  # noqa: C901
    database: REQ_ARG_TYPE_DATABASE,
    outdir: REQ_ARG_TYPE_OUTDIR,
    run_id: OPT_ARG_TYPE_RUN_ID = None,
    # Would like to replace this with Literal["md5", "filename", "stem"] once typer updated
    label: Annotated[
        str,
        typer.Option(
            click_type=click.Choice(["md5", "filename", "stem"]),
            help="How to label the genomes",
        ),
    ] = "stem",
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
        print(f"INFO: Plotting run-id {run_id}")

    method = run.configuration.method

    heatmaps_done = 0
    for matrix, name in (
        (run.identities, "identity"),
        (run.cov_query, "query_cov"),
        (run.hadamard, "hadamard"),
    ):
        if matrix is None:
            # This is mainly for mypy to assert the matrix is not None
            msg = f"ERROR: Could not load run {method} matrix"  # pragma: no cover
            sys.exit(msg)  # pragma: no cover

        nulls = int(matrix.isnull().sum().sum())
        n = len(matrix)
        if nulls:
            msg = (
                f"WARNING: Cannot plot {name} as matrix contains {nulls} nulls"
                f" (out of {n}²={n**2} {method} comparisons)\n"
            )
            sys.stderr.write(msg)
            continue

        try:
            matrix = run.relabelled_matrix(matrix, label)  # noqa: PLW2901
        except ValueError as err:
            msg = f"ERROR: {err}"
            sys.exit(msg)

        # Next cluster the matrix & prepare the figure (move code to new file)
        import warnings

        import seaborn as sns

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=(
                    "scipy.cluster: The symmetric non-negative hollow observation"
                    " matrix looks suspiciously like an uncondensed distance matrix"
                ),
            )
            figure = sns.clustermap(matrix)

        for ext in ("tsv", "png", "jpg", "svg", "pdf"):
            filename = outdir / f"{method}_{name}.{ext}"
            if ext == "tsv":
                # Apply the clustering reordering to match the figure:
                matrix = matrix.iloc[  # noqa: PLW2901
                    figure.dendrogram_row.reordered_ind,
                    figure.dendrogram_row.reordered_ind,
                ]
                matrix.to_csv(filename, sep="\t")
            else:
                figure.savefig(filename)

        heatmaps_done += 1
        # Next want to plot distributions of the scores (scatter plots)

    if not heatmaps_done:
        msg = "ERROR: Unable to plot any heatmaps (check for nulls)"
        sys.exit(msg)
    print(f"Wrote {heatmaps_done} heatmaps to {outdir}/{method}_*.*")

    session.close()
    return 0


if __name__ == "__main__":
    sys.exit(app())  # pragma: no cover
