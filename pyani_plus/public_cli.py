# The MIT License
#
# Copyright (c) 2024 University of Strathclyde
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

import multiprocessing
import sys
import tempfile
from pathlib import Path
from time import sleep
from typing import Annotated

import pandas as pd
import typer
from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table
from rich.text import Text
from sqlalchemy import insert, text
from sqlalchemy.exc import NoResultFound, OperationalError
from sqlalchemy.orm import Session

from pyani_plus import FASTA_EXTENSIONS, db_orm, tools
from pyani_plus.methods import method_anib, method_fastani
from pyani_plus.utils import available_cores, file_md5sum
from pyani_plus.workflows import SnakemakeRunner

# Reused required command line arguments (which have no default)
REQ_ARG_TYPE_DATABASE = Annotated[
    Path,
    typer.Option(
        help="Path to pyANI-plus SQLite3 database",
        show_default=False,
        dir_okay=False,
        file_okay=True,
    ),
]
REQ_ARG_TYPE_RUN_NAME = Annotated[
    str, typer.Option(help="Run name", show_default=False)
]
REQ_ARG_TYPE_OUTDIR = Annotated[
    Path,
    typer.Option(
        help="Output directory",
        show_default=False,
        exists=True,
        dir_okay=True,
        file_okay=False,
    ),
]
REQ_ARG_TYPE_FASTA_DIR = Annotated[
    Path,
    typer.Argument(
        help=f"Directory of FASTA files (extensions {', '.join(sorted(FASTA_EXTENSIONS))})",
        show_default=False,
    ),
]

# Reused optional command line arguments (defined with a default):
OPT_ARG_TYPE_FRAGSIZE = Annotated[
    int,
    typer.Option(
        help="Comparison method fragment size",
        rich_help_panel="Method parameters",
        min=1,
    ),
]
# fastANI has maximum (and default) k-mer size 16
OPT_ARG_TYPE_KMERSIZE = Annotated[
    int,
    typer.Option(
        help="Comparison method k-mer size", rich_help_panel="Method parameters", min=1
    ),
]
OPT_ARG_TYPE_MINMATCH = Annotated[
    float,
    typer.Option(
        help="Comparison method min-match",
        rich_help_panel="Method parameters",
        min=0.0,
        max=1.0,
    ),
]
OPT_ARG_TYPE_CREATE_DB = Annotated[
    # Listing name(s) explicitly to avoid automatic matching --no-create-db
    bool, typer.Option("--create-db", help="Create database if does not exist")
]

progress_columns = [
    TextColumn("[progress.description]{task.description}"),
    BarColumn(),
    TaskProgressColumn(),
    # Removing TimeRemainingColumn() from defaults, replacing with:
    TimeElapsedColumn(),
    # Add this last as have some out of N and some out of N^2:
    MofNCompleteColumn(),
]

app = typer.Typer(
    context_settings={"help_option_names": ["-h", "--help"]},
)


def check_db(database: Path | str, create_db: bool) -> None:  # noqa: FBT001
    """Check DB exists, or using create_db=True."""
    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist, but not using --create-db"
        sys.exit(msg)


def check_fasta(fasta: Path) -> list[Path]:
    """Check fasta is a directory and return list of FASTA files in it."""
    if not fasta.is_dir():
        msg = f"ERROR: FASTA input {fasta} is not a directory"
        sys.exit(msg)

    fasta_names: list[Path] = []
    for pattern in FASTA_EXTENSIONS:
        fasta_names.extend(fasta.glob("*" + pattern))
    if not fasta_names:
        msg = f"ERROR: No FASTA input genomes under {fasta} with extensions {', '.join(FASTA_EXTENSIONS)}"
        sys.exit(msg)

    return fasta_names


def progress_bar_via_db_comparisons(
    database: Path, run_id: int, interval: float = 1.0
) -> None:
    """Show a progress bar based on monitoring the DB entries.

    Will only self-terminate once all the run's comparisons are
    in the database! Up to the caller to terminate it early.
    """
    session = db_orm.connect_to_db(database)
    run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()

    already_done = run.comparisons().count()
    total = run.genomes.count() ** 2 - already_done

    done = 0
    old_db_version = -1
    with Progress(*progress_columns) as progress:
        task = progress.add_task("Comparing pairs", total=total)
        while done < total:
            sleep(interval)
            # Have there been any DB changes?
            try:
                db_version = (
                    session.connection().execute(text("PRAGMA data_version;")).one()[0]
                )
            except OperationalError:
                # Probably another process is updating the DB, try again soon
                continue
            if old_db_version == db_version:
                # No changes, check again soon
                continue
            old_db_version = db_version
            new = run.comparisons().count() - already_done - done
            progress.update(task, advance=new)
            done += new
    session.close()


def run_snakemake_with_progress_bar(  # noqa: PLR0913
    workflow_name: str,
    database: Path,
    run_id: int,
    targets: list[Path],
    params: dict,
    working_directory: Path,
    interval: float = 0.5,
) -> None:
    """Run snakemake with a progress bar."""
    runner = SnakemakeRunner(workflow_name)
    # All rest replaces one line, runner.run_workflow(targets, params, workdir=working_directory)

    # As of Python 3.8 onwards, the default on macOS ("Darwin") is "spawn"
    # As of Python 3.12, the default of "fork" on Linux triggers a deprecation warning.
    # This should match the defaults on Python 3.14 onwards.
    # Note mypy currently can't handle this dynamic situation, their issue #8603
    p = multiprocessing.get_context(  # type:ignore [attr-defined]
        "spawn" if sys.platform == "darwin" else "forkserver"
    ).Process(
        target=progress_bar_via_db_comparisons,
        args=(
            database,
            run_id,
            interval,
        ),
    )
    p.start()

    # Call slurm!
    try:
        runner.run_workflow(targets, params, workdir=working_directory)
    except Exception as err:  # noqa: BLE001
        p.terminate()
        # Ensure traceback starts on a new line after interrupted progress bar
        print()  # noqa: T201
        raise err from None

    if p.is_alive():
        # Should have finished, but perhaps final update is pending...
        sleep(interval)
        p.terminate()
    p.join()


def record_genomes(
    session: Session, run: db_orm.Run, fasta_names: list[Path]
) -> dict[Path, str]:
    """Compute the MD5 of the FASTA files and record them against this run.

    Returns a dict mapping filenames to MD5 checksums.
    """
    # Reuse existing genome entries and/or log new ones
    n = len(fasta_names)
    filename_to_md5 = {}
    hashes = set()
    with Progress(*progress_columns) as progress:
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

    # Will now update the runs_genomes linker table, but did not scale
    # well via runs.genomes = [list of genome ORM objects]
    # We could update the runs_genomes incrementally, or in batches?
    run_id = run.run_id
    session.execute(
        insert(db_orm.RunGenomeAssociation),
        [{"run_id": run_id, "genome_hash": md5} for md5 in hashes],
    )
    del hashes
    run.status = "Setup"
    session.commit()
    print(f"Run setup with {n} genomes in database")  # noqa: T201
    return filename_to_md5


def run_method(  # noqa: PLR0913
    database: Path,
    name: str,
    method: str,
    fasta: Path,
    target_extension: str,
    tool: tools.ExternalToolData,
    fragsize: int | None,
    maxmatch: bool | None,
    kmersize: int | None,
    minmatch: float | None,
    params: dict[str, object],
) -> int:
    """Run the snakemake workflow for given method and log run to database."""
    fasta_names = check_fasta(fasta)

    workflow_name = f"snakemake_{method.lower()}.smk"

    # We should have caught all the obvious failures above,
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
        maxmatch,
        kmersize,
        minmatch,
        create=True,
    )

    run = db_orm.add_run(
        session,
        config,
        cmdline=" ".join(sys.argv),
        status="Initialising",
        name=name,
        date=None,
        genomes=[],
    )
    run_id = run.run_id

    n = len(fasta_names)
    filename_to_md5 = record_genomes(session, run, fasta_names)
    del fasta_names

    done = run.comparisons().count()
    if done == n**2:
        print(  # noqa: T201
            f"Database already has all {n}x{n}={n**2} comparisons"
        )
    else:
        print(  # noqa: T201
            f"Database already has {done} of {n}x{n}={n**2} comparisons, {n**2 - done} needed"
        )
        done_hashes = {(_.query_hash, _.subject_hash) for _ in run.comparisons()}
        session.close()  # Reduce chance of DB locking

        # Run snakemake wrapper
        with tempfile.TemporaryDirectory(prefix="pyani-plus_") as tmp:
            tmp_path = Path(tmp) / "working"
            out_path = Path(tmp) / "output"
            targets = [
                out_path / f"{query.stem}_vs_{subject.stem}{target_extension}"
                for query in filename_to_md5
                for subject in filename_to_md5
                if (filename_to_md5[query], filename_to_md5[subject]) not in done_hashes
            ]
            del done_hashes
            # Must all inputs be in one place? Symlinks?
            params["indir"] = fasta.resolve()  # must be absolute
            params["outdir"] = out_path.resolve()
            params["db"] = Path(database).resolve()  # must be absolute
            params["cores"] = available_cores()  # should make configurable
            run_snakemake_with_progress_bar(
                workflow_name,
                Path(database),
                run_id,
                targets,
                params,
                tmp_path,
            )

        # Reconnect to the DB
        session = db_orm.connect_to_db(database)
        run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()
        done = run.comparisons().count()

    if done != n**2:
        msg = f"ERROR: Only have {done} of {n}x{n}={n**2} comparisons needed"
        sys.exit(msg)

    run.cache_comparisons()  # will this needs a progress bar too with large n?
    run.status = "Done"
    session.commit()
    session.close()
    print(  # noqa: T201
        f"Logged new run-id {run_id} for method {method} with {n} genomes to database {database}"
    )
    return 0


@app.command(rich_help_panel="ANI methods")
def anim(
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the run table:
    name: REQ_ARG_TYPE_RUN_NAME,
    # Does not use fragsize, maxmatch, kmersize, or minmatch
    create_db: OPT_ARG_TYPE_CREATE_DB = False,  # noqa: FBT002
) -> int:
    """Execute ANIm calculations, logged to a pyANI-plus SQLite3 database."""
    check_db(database, create_db)

    target_extension = ".filter"
    tool = tools.get_nucmer()
    params = {
        "nucmer": tool.exe_path,
        "delta_filter": tools.get_delta_filter().exe_path,
        "mode": "mum",
    }
    fragsize = maxmatch = kmersize = minmatch = None
    return run_method(
        database,
        name,
        "ANIm",
        fasta,
        target_extension,
        tool,
        fragsize,
        maxmatch,
        kmersize,
        minmatch,
        params,
    )


@app.command(rich_help_panel="ANI methods")
def dnadiff(
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the run table:
    name: REQ_ARG_TYPE_RUN_NAME,
    # Does not use fragsize, maxmatch, kmersize, or minmatch
    create_db: OPT_ARG_TYPE_CREATE_DB = False,  # noqa: FBT002
) -> int:
    """Execute mumer-based dnadiff calculations, logged to a pyANI-plus SQLite3 database."""
    check_db(database, create_db)

    target_extension = ".qdiff"  # or .mcoords as rule makes both
    # We don't actually call the tool dnadiff (which has its own version),
    # rather we call nucmer, delta-filter, show-diff and show-coords from MUMmer
    tool = tools.get_nucmer()
    params: dict[str, object] = {
        "nucmer": tool.exe_path,
        "delta_filter": tools.get_delta_filter().exe_path,
        "show_diff": tools.get_show_diff().exe_path,
        "show_coords": tools.get_show_coords().exe_path,
    }
    fragsize = maxmatch = kmersize = minmatch = None
    return run_method(
        database,
        name,
        "dnadiff",
        fasta,
        target_extension,
        tool,
        fragsize,
        maxmatch,
        kmersize,
        minmatch,
        params,
    )


@app.command(rich_help_panel="ANI methods")
def anib(
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the run table:
    name: REQ_ARG_TYPE_RUN_NAME,
    # These are all for the configuration table:
    fragsize: OPT_ARG_TYPE_FRAGSIZE = method_anib.FRAGSIZE,
    # Does not use maxmatch, kmersize, or minmatch
    create_db: OPT_ARG_TYPE_CREATE_DB = False,  # noqa: FBT002
) -> int:
    """Execute ANIb calculations, logged to a pyANI-plus SQLite3 database."""
    check_db(database, create_db)

    target_extension = ".tsv"
    tool = tools.get_blastn()
    alt = tools.get_makeblastdb()
    if tool.version != alt.version:
        msg = f"ERROR: blastn {tool.version} vs makeblastdb {alt.version}"
        sys.exit(msg)
    params = {
        "blastn": tool.exe_path,
        "makeblastdb": alt.exe_path,
        "fragLen": fragsize,
    }
    maxmatch = kmersize = minmatch = None
    return run_method(
        database,
        name,
        "ANIb",
        fasta,
        target_extension,
        tool,
        fragsize,
        maxmatch,
        kmersize,
        minmatch,
        params,
    )


@app.command(rich_help_panel="ANI methods")
def fastani(  # noqa: PLR0913
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the run table:
    name: REQ_ARG_TYPE_RUN_NAME,
    # These are all for the configuration table:
    fragsize: OPT_ARG_TYPE_FRAGSIZE = method_fastani.FRAG_LEN,
    # Does not use maxmatch
    kmersize: OPT_ARG_TYPE_KMERSIZE = method_fastani.KMER_SIZE,
    minmatch: OPT_ARG_TYPE_MINMATCH = method_fastani.MIN_FRACTION,
    create_db: OPT_ARG_TYPE_CREATE_DB = False,  # noqa: FBT002
) -> int:
    """Execute fastANI calculations, logged to a pyANI-plus SQLite3 database."""
    check_db(database, create_db)

    target_extension = ".fastani"
    tool = tools.get_fastani()
    params = {
        "fastani": tool.exe_path,
        "fragLen": fragsize,
        "kmerSize": kmersize,
        "minFrac": minmatch,
    }
    maxmatch = None
    return run_method(
        database,
        name,
        "fastANI",
        fasta,
        target_extension,
        tool,
        fragsize,
        maxmatch,
        kmersize,
        minmatch,
        params,
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
    table.add_column(Text("Comparisons", justify="left"), justify="right", no_wrap=True)
    table.add_column("Status")
    table.add_column("Name")
    # Would be nice to report {conf.program} {conf.version} too,
    # perhaps conditional on the terminal width?
    for run in runs:
        conf = run.configuration
        n = run.genomes.count()
        done = run.comparisons().count()
        comparison_summary = f"{done}/{n**2}={n}²"
        if done == 0:
            comparison_summary = Text(comparison_summary, style="red")
        elif done != n**2:
            comparison_summary = Text(comparison_summary, style="yellow")
        else:
            comparison_summary = Text(comparison_summary, style="green")
        table.add_row(
            str(run.run_id),
            str(run.date.date()),
            conf.method,
            comparison_summary,
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
    run_id: Annotated[
        int | None,
        typer.Option(help="Which run to report (optional if DB contains only one)"),
    ] = None,
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

    if run_id is None:
        runs = session.query(db_orm.Run)
        count = runs.count()
        if count == 1:
            run = runs.one()
            run_id = run.run_id
            print(f"INFO: Reporting on run-id {run_id} from {database}")  # noqa: T201
        elif count:
            msg = (
                f"ERROR: Database {database} contains {count} runs,"
                " use --run-id to specify which."
                " Use the list-runs command for more information."
            )
            sys.exit(msg)
        else:
            msg = f"ERROR: Database {database} contains no runs."
            sys.exit(msg)
    else:
        try:
            run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()
        except NoResultFound:
            msg = (
                f"ERROR: Database {database} has no run-id {run_id}."
                " Use the list-runs command for more information."
            )
            sys.exit(msg)

    conf = run.configuration

    if not run.comparisons().count():
        msg = f"ERROR: Database {database} run-id {run_id} has no comparisons"
        sys.exit(msg)
    # What if the run is incomplete? Just output with NaN?
    if run.identities is None:
        run.cache_comparisons()
    if not isinstance(run.identities, pd.DataFrame):
        msg = f"ERROR: Could not access identities matrix from JSON for {run_id}"
        sys.exit(msg)

    # Question: Should we export a plain text of JSON summary of the configuration etc?
    # Question: Should we include the run-id in the filenames name?
    # Question: Should we match the property in filenames to old pyANI (e.g. coverage)?
    # Question: Should we offer MD5 alternatives, and how? e.g. filename or labels
    # (seems more suited to a report command offering tables and plots)
    run.identities.to_csv(outdir / f"{conf.method}_identity.tsv", sep="\t")
    run.aln_length.to_csv(outdir / f"{conf.method}_aln_lengths.tsv", sep="\t")
    run.sim_errors.to_csv(outdir / f"{conf.method}_sim_errors.tsv", sep="\t")
    run.cov_query.to_csv(outdir / f"{conf.method}_query_cov.tsv", sep="\t")
    run.hadamard.to_csv(outdir / f"{conf.method}_hadamard.tsv", sep="\t")

    print(f"Wrote matrices to {outdir}/{conf.method}_*.tsv")  # noqa: T201

    session.close()
    return 0


if __name__ == "__main__":
    sys.exit(app())  # pragma: no cover
