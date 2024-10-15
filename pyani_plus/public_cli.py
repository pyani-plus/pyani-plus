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

import sys
import tempfile
from pathlib import Path
from typing import Annotated

import pandas as pd
import typer
from rich.progress import track
from sqlalchemy.exc import NoResultFound

from pyani_plus import db_orm, tools
from pyani_plus.snakemake import snakemake_scheduler
from pyani_plus.utils import file_md5sum

FASTA_EXTENSIONS = {".fasta", ".fas", ".fna"}  # define more centrally?

app = typer.Typer()


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

    # Reuse existing genome entries and/or log new ones
    genomes = []
    filename_to_md5 = {}
    for filename in track(
        fasta_names, description=f"Processing {len(fasta_names)} FASTA"
    ):
        md5 = file_md5sum(filename)
        genomes.append(db_orm.db_genome(session, filename, md5, create=True))
        filename_to_md5[filename] = md5
    n = len(genomes)

    run = db_orm.add_run(
        session,
        config,
        cmdline=" ".join(sys.argv),
        status="Setup",
        name=name,
        date=None,
        genomes=genomes,
    )
    run_id = run.run_id
    session.commit()  # Redundant?
    del genomes

    done = run.comparisons().count()
    if done < n**2:
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
            params["cores"] = 1  # should make dynamic
            runner = snakemake_scheduler.SnakemakeRunner(workflow_name)
            runner.run_workflow(targets, params, workdir=tmp_path)

        # Reconnect to the DB
        session = db_orm.connect_to_db(database)
        run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()
        done = run.comparisons().count()

    if done != n**2:
        msg = f"ERROR: Only have {done} of {n}x{n}={n**2} comparisons needed"
        sys.exit(msg)

    run.cache_comparisons()
    run.status = "Done"
    session.commit()
    session.close()
    print(  # noqa: T201
        f"Logged new run-id {run_id} for method {method} with {n} genomes to database {database}"
    )
    return 0


@app.command()
def anim(
    fasta: Annotated[
        Path,
        typer.Argument(
            help=f"Directory of FASTA files (extensions {', '.join(sorted(FASTA_EXTENSIONS))})",
            show_default=False,
        ),
    ],
    database: Annotated[
        Path,
        typer.Option(help="Path to pyANI-plus SQLite3 database", show_default=False),
    ],
    # These are for the run table:
    name: Annotated[str, typer.Option(help="Run name", show_default=False)],
    # Does not use fragsize, maxmatch, kmersize, or minmatch
    create_db: Annotated[  # noqa: FBT002
        bool, typer.Option(help="Create database if does not exist")
    ] = False,
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


@app.command()
def anib(
    fasta: Annotated[
        Path,
        typer.Argument(
            help=f"Directory of FASTA files (extensions {', '.join(sorted(FASTA_EXTENSIONS))})",
            show_default=False,
        ),
    ],
    database: Annotated[
        Path,
        typer.Option(help="Path to pyANI-plus SQLite3 database", show_default=False),
    ],
    # These are for the run table:
    name: Annotated[str, typer.Option(help="Run name", show_default=False)],
    # These are all for the configuration table:
    fragsize: Annotated[
        int, typer.Option(help="Comparison method fragment size")
    ] = 1020,
    # Does not use maxmatch, kmersize, or minmatch
    create_db: Annotated[  # noqa: FBT002
        bool, typer.Option(help="Create database if does not exist")
    ] = False,
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


@app.command()
def fastani(  # noqa: PLR0913
    fasta: Annotated[
        Path,
        typer.Argument(
            help=f"Directory of FASTA files (extensions {', '.join(sorted(FASTA_EXTENSIONS))})",
            show_default=False,
        ),
    ],
    database: Annotated[
        Path,
        typer.Option(help="Path to pyANI-plus SQLite3 database", show_default=False),
    ],
    # These are for the run table:
    name: Annotated[str, typer.Option(help="Run name", show_default=False)],
    # These are all for the configuration table:
    fragsize: Annotated[
        int, typer.Option(help="Comparison method fragment size")
    ] = 3000,
    # Does not use maxmatch
    kmersize: Annotated[int, typer.Option(help="Comparison method k-mer size")] = 16,
    minmatch: Annotated[float, typer.Option(help="Comparison method min-match")] = 0.2,
    create_db: Annotated[  # noqa: FBT002
        bool, typer.Option(help="Create database if does not exist")
    ] = False,
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
    database: Annotated[str, typer.Option(help="Path to pyANI-plus SQLite3 database")],
) -> int:
    """List the runs defined in a given pyANI-plus SQLite3 database."""
    if database == ":memory:" or not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    session = db_orm.connect_to_db(database)
    runs = session.query(db_orm.Run)
    print(f"Database {database} contains {runs.count()} runs")  # noqa: T201
    for run in runs:
        print(f"Run-ID {run.run_id}, '{run.name}'")  # noqa: T201
        print(f"  {run.genomes.count()} genomes; date {run.date}; status {run.status}")  # noqa: T201
        conf = run.configuration
        print(f"  Method {conf.method}; program {conf.program} {conf.version}")  # noqa: T201
    session.close()
    return 0


@app.command()
def export_run(
    database: Annotated[str, typer.Option(help="Path to pyANI-plus SQLite3 database")],
    outdir: Annotated[Path, typer.Option(help="Output directory")],
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
