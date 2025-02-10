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
"""Implements the private command line interface (CLI) used internally.

The commands defined here are intended to be used from within pyANI-plus via
snakemake, for example from worker nodes, to log results to the database.
"""

import os
import platform
import signal
import sys
import tempfile
from contextlib import nullcontext
from pathlib import Path
from typing import Annotated

import typer
from rich.progress import Progress
from sqlalchemy.dialects.sqlite import insert as sqlite_insert
from sqlalchemy.orm import Session

from pyani_plus import PROGRESS_BAR_COLUMNS, db_orm, tools
from pyani_plus.methods import anib, anim, dnadiff, fastani, sourmash
from pyani_plus.public_cli_args import (
    OPT_ARG_TYPE_CREATE_DB,
    OPT_ARG_TYPE_TEMP,
    REQ_ARG_TYPE_DATABASE,
    REQ_ARG_TYPE_FASTA_DIR,
)
from pyani_plus.utils import (
    check_fasta,
    check_output,
    file_md5sum,
    filename_stem,
    stage_file,
)

app = typer.Typer(
    context_settings={"help_option_names": ["-h", "--help"]},
)

REQ_ARG_TYPE_RUN_ID = Annotated[
    int,
    typer.Option(help="Database run ID", show_default=False),
]
REQ_ARG_TYPE_CONFIG_ID = Annotated[
    int,
    typer.Option(help="Database configuration ID", show_default=False),
]
OPT_ARG_TYPE_QUIET = Annotated[
    # Listing name(s) explicitly to avoid automatic matching --no-quiet
    bool, typer.Option("--quiet", help="Suppress any output except if fails")
]
REQ_ARG_TYPE_RUN_NAME = Annotated[
    str, typer.Option(help="Run name", show_default=False)
]
REQ_ARG_TYPE_FASTA_FILES = Annotated[
    list[Path],
    typer.Argument(
        help="Path(s) to FASTA file(s)",
        show_default=False,
        exists=True,
        dir_okay=False,
        file_okay=True,
    ),
]
REQ_ARG_TYPE_QUERY_FASTA = Annotated[
    Path,
    typer.Option(
        help="Path to query FASTA file",
        show_default=False,
        exists=True,
        dir_okay=False,
        file_okay=True,
    ),
]
REQ_ARG_TYPE_SUBJECT_FASTA = Annotated[
    Path,
    typer.Option(
        help="Path to subject (reference) FASTA file",
        show_default=False,
        exists=True,
        dir_okay=False,
        file_okay=True,
    ),
]
REQ_ARG_TYPE_SUBJECT = Annotated[
    str,
    typer.Option(
        help="Filename or hash of subject (reference) FASTA file",
        show_default=False,
    ),
]
REQ_ARG_TYPE_METHOD = Annotated[
    str, typer.Option(help="Method, e.g. ANIm", show_default=False)
]
REQ_ARG_TYPE_PROGRAM = Annotated[
    str, typer.Option(help="Program, e.g. nucmer", show_default=False)
]
REQ_ARG_TYPE_VERSION = Annotated[
    str, typer.Option(help="Program version, e.g. 3.1", show_default=False)
]

# Reused optional command line arguments (used here with None as their default,
# whereas in the public_cli they have type-appropriate defaults):
NONE_ARG_TYPE_FRAGSIZE = Annotated[
    int | None,
    typer.Option(
        help="Comparison method fragment size",
        rich_help_panel="Method parameters",
        min=1,
    ),
]
NONE_ARG_TYPE_MODE = Annotated[
    str | None,
    typer.Option(
        help="Comparison method specific mode", rich_help_panel="Method parameters"
    ),
]
NONE_ARG_TYPE_KMERSIZE = Annotated[
    int | None,
    typer.Option(
        help="Comparison method k-mer size", rich_help_panel="Method parameters", min=1
    ),
]
NONE_ARG_TYPE_MINMATCH = Annotated[
    float | None,
    typer.Option(
        help="Comparison method min-match",
        rich_help_panel="Method parameters",
        min=0.0,
        max=1.0,
    ),
]
NONE_ARG_TYPE_EXTRA = Annotated[
    str | None,
    typer.Option(
        help="Comparison method specific mode", rich_help_panel="Method parameters"
    ),
]


def _check_tool_version(
    tool: tools.ExternalToolData, configuration: db_orm.Configuration
) -> None:
    """Confirm the tool and version matches the given configuration."""
    if (
        configuration.program != tool.exe_path.stem
        or configuration.version != tool.version
    ):
        msg = (
            "ERROR: Run configuration was"
            f" {configuration.program} {configuration.version}"
            f" but we have {tool.exe_path.stem} {tool.version}"
        )
        sys.exit(msg)


@app.command(rich_help_panel="Low-level logging")
def log_configuration(  # noqa: PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    method: REQ_ARG_TYPE_METHOD,
    program: REQ_ARG_TYPE_PROGRAM,
    version: REQ_ARG_TYPE_VERSION,
    *,
    fragsize: NONE_ARG_TYPE_FRAGSIZE = None,
    mode: NONE_ARG_TYPE_MODE = None,
    kmersize: NONE_ARG_TYPE_KMERSIZE = None,
    minmatch: NONE_ARG_TYPE_MINMATCH = None,
    extra: NONE_ARG_TYPE_EXTRA = None,
    create_db: OPT_ARG_TYPE_CREATE_DB = False,
) -> int:
    """Log a specific method configuration to database.

    Any pre-existing configuration entry is left as is.
    """
    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist, but not using --create-db"
        sys.exit(msg)

    print(f"Logging configuration to {database}")
    session = db_orm.connect_to_db(database)
    config = db_orm.db_configuration(
        session=session,
        method=method,
        program=program,
        version=version,
        fragsize=fragsize,
        mode=mode,
        kmersize=kmersize,
        minmatch=minmatch,
        extra=extra,
        create=True,
    )
    session.commit()  # should be redundant
    print(f"Configuration identifier {config.configuration_id}")
    session.close()

    return 0


@app.command(rich_help_panel="Low-level logging")
def log_genome(
    fasta: REQ_ARG_TYPE_FASTA_FILES,
    database: REQ_ARG_TYPE_DATABASE,
    *,
    create_db: OPT_ARG_TYPE_CREATE_DB = False,
) -> int:
    """Compute MD5 checksums of given FASTA files, log them to database.

    Any pre-existing duplicate FASTA entries are left as is.
    """
    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist, but not using --create-db"
        sys.exit(msg)

    print(f"Logging genome to {database}")
    session = db_orm.connect_to_db(database)

    file_total = 0
    if fasta:
        with Progress(*PROGRESS_BAR_COLUMNS) as progress:
            for filename in progress.track(fasta, description="Processing..."):
                file_total += 1
                md5 = file_md5sum(filename)
                db_orm.db_genome(session, filename, md5, create=True)
    session.commit()
    session.close()
    print(f"Processed {file_total} FASTA files")

    return 0


@app.command(rich_help_panel="Low-level logging")
def log_run(  # noqa: PLR0913
    fasta: REQ_ARG_TYPE_FASTA_DIR,
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the run table:
    cmdline: Annotated[str, typer.Option(help="Run command line", show_default=False)],
    status: Annotated[str, typer.Option(help="Run status", show_default=False)],
    name: REQ_ARG_TYPE_RUN_NAME,
    # These are all for the configuration table:
    method: REQ_ARG_TYPE_METHOD,
    program: REQ_ARG_TYPE_PROGRAM,
    version: REQ_ARG_TYPE_VERSION,
    *,
    fragsize: NONE_ARG_TYPE_FRAGSIZE = None,
    mode: NONE_ARG_TYPE_MODE = None,
    kmersize: NONE_ARG_TYPE_KMERSIZE = None,
    extra: NONE_ARG_TYPE_EXTRA = None,
    minmatch: NONE_ARG_TYPE_MINMATCH = None,
    create_db: OPT_ARG_TYPE_CREATE_DB = False,
) -> int:
    """Log a run (and if need be, associated configuration and genome rows).

    There is currently no easy way to update an existing run (e.g. once more
    comparisons have been completed and you want to refresh the cached matrices
    and update the run status).
    """
    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist, but not using --create-db"
        sys.exit(msg)

    print(f"Logging run to {database}")
    session = db_orm.connect_to_db(database)

    # Reuse existing config, or log a new one
    config = db_orm.db_configuration(
        session=session,
        method=method,
        program=program,
        version=version,
        fragsize=fragsize,
        mode=mode,
        kmersize=kmersize,
        minmatch=minmatch,
        extra=extra,
        create=True,
    )

    fasta_to_hash = {}
    fasta_names = check_fasta(fasta)
    if fasta_names:
        # Reuse existing genome entries and/or log new ones
        with Progress(*PROGRESS_BAR_COLUMNS) as progress:
            for filename in progress.track(fasta_names, description="Processing..."):
                md5 = file_md5sum(filename)
                fasta_to_hash[filename] = md5
                db_orm.db_genome(session, filename, md5, create=True)

    run = db_orm.add_run(
        session,
        config,
        cmdline,
        fasta,
        status,
        name,
        date=None,
        fasta_to_hash=fasta_to_hash,
    )
    # No point caching empty matrices, even partial ones is debatable
    if run.comparisons().count() == len(fasta_to_hash) ** 2:
        run.cache_comparisons()
    run_id = run.run_id

    session.commit()
    session.close()
    print(f"Run identifier {run_id}")

    return 0


@app.command(rich_help_panel="Low-level logging")
def log_comparison(  # noqa: PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
    config_id: REQ_ARG_TYPE_CONFIG_ID,
    query_fasta: REQ_ARG_TYPE_QUERY_FASTA,
    subject_fasta: REQ_ARG_TYPE_SUBJECT_FASTA,
    *,
    # Optional comparison table entries
    aln_length: Annotated[int | None, typer.Option(help="Alignment length")] = None,
    identity: Annotated[float | None, typer.Option(help="Percent identity")] = None,
    sim_errors: Annotated[int | None, typer.Option(help="Alignment length")] = None,
    cov_query: Annotated[float | None, typer.Option(help="Alignment length")] = None,
    cov_subject: Annotated[float | None, typer.Option(help="Alignment length")] = None,
) -> int:
    """Log single pairwise comparison to database."""
    if database != ":memory:" and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    print(f"Logging comparison to {database}")
    session = db_orm.connect_to_db(database)
    # Give a better error message that if adding comparison fails:
    if (
        not session.query(db_orm.Configuration)
        .where(db_orm.Configuration.configuration_id == config_id)
        .count()
    ):
        msg = f"ERROR - {database} does not contain configuration_id={config_id}"
        sys.exit(msg)

    query_md5 = file_md5sum(query_fasta)
    db_orm.db_genome(session, query_fasta, query_md5)

    subject_md5 = file_md5sum(subject_fasta)
    db_orm.db_genome(session, subject_fasta, subject_md5)

    db_orm.db_comparison(
        session,
        configuration_id=config_id,
        query_hash=query_md5,
        subject_hash=subject_md5,
        identity=identity,
        aln_length=aln_length,
        sim_errors=sim_errors,
        cov_query=cov_query,
        cov_subject=cov_subject,
    )

    session.commit()
    return 0


@app.command()
def compute_column(  # noqa: C901
    database: REQ_ARG_TYPE_DATABASE,
    run_id: REQ_ARG_TYPE_RUN_ID,
    subject: Annotated[
        str,
        typer.Option(
            help="Subject (reference) FASTA filename, MD5 checksum, or index (integer).",
            show_default=False,
            exists=True,
            dir_okay=False,
            file_okay=True,
        ),
    ],
    *,
    temp: OPT_ARG_TYPE_TEMP = None,
    quiet: OPT_ARG_TYPE_QUIET = False,
) -> int:
    """Run the method for one column and log pairwise comparisons to the database.

    The database and run ID specify the method and configuration, additionally
    you must supply a subject filename, hash, or column index to control which
    column of the matrix is to be computed.

    If using a column number, these are taken to be zero based but it will accept
    0 or N to mean the first subject. This is intended to facilitate use with
    cluster array jobs.
    """
    if database != ":memory:" and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    # We want to receive any SIGINT as a KeyboardInterrupt even if we
    # are run via a bash shell or other non-interactive setting:
    signal.signal(signal.SIGINT, signal.default_int_handler)

    # Under SLURM when a job is cancelled via scancel, it sends SIGTERM
    # for a graceful shutdown then 30s later a SIGKILL hard kill.
    # For simplicity, treat SIGTERM as a KeyboardInterrupt too.
    signal.signal(signal.SIGTERM, signal.default_int_handler)

    session = db_orm.connect_to_db(database)
    run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()
    config = run.configuration
    method = config.method

    filename_to_hash = {_.fasta_filename: _.genome_hash for _ in run.fasta_hashes}
    hash_to_filename = {_.genome_hash: _.fasta_filename for _ in run.fasta_hashes}
    n = len(hash_to_filename)

    if subject in hash_to_filename:
        subject_hash = subject
    elif Path(subject).name in filename_to_hash:
        subject_hash = filename_to_hash[Path(subject).name]
    else:
        try:
            column = int(subject)
        except ValueError:
            msg = f"ERROR: Did not recognise {subject!r} as an MD5 hash, filename, or column number in run-id {run_id}"
            sys.exit(msg)
        if column < 0 or len(hash_to_filename) < column:
            msg = f"ERROR: Column should be in range 0 to {n}, not {subject}"
            sys.exit(msg)
        if column == n:
            if not quiet:
                sys.stderr.write("INFO: Treating subject N as 0 (first column)\n")
            column = 0
        subject_hash = sorted(hash_to_filename)[column]
        del column

    # What comparisons are needed?
    query_hashes: list[str] = sorted(
        set(hash_to_filename).difference(
            comp.query_hash
            for comp in run.comparisons().where(
                db_orm.Comparison.subject_hash == subject_hash
            )
        )
    )

    if not query_hashes:
        if not quiet:
            print(f"INFO: No {method} comparisons needed against {subject_hash}")
        return 0

    # Will probably want to move each of these functions to the relevant method module...
    try:
        compute = {
            "fastANI": compute_fastani,
            "ANIb": compute_anib,
            "ANIm": compute_anim,
            "dnadiff": compute_dnadiff,
            "external-alignment": compute_external_alignment,
        }[method]
    except KeyError:
        msg = f"ERROR: Unknown method {method} for run-id {run_id} in {database}"
        sys.exit(msg)

    # On a cluster we are likely in a temp working directory, meaning
    # if it is a relative path, run.fasta_directory is useless without
    # evaluating it relative to the DB filename.
    fasta_dir = Path(run.fasta_directory)
    if not fasta_dir.is_absolute():
        fasta_dir = (database.parent / fasta_dir).absolute()

    # Either use the specified temp-directory (and do not clean up),
    # or use a system temp-directory (and do clean up)
    with nullcontext(temp) if temp else tempfile.TemporaryDirectory() as tmp_dir:
        return compute(
            Path(tmp_dir),
            session,
            run,
            fasta_dir,
            hash_to_filename,
            filename_to_hash,
            query_hashes,
            subject_hash,
            quiet=quiet,
        )


def compute_fastani(  # noqa: PLR0913
    tmp_dir: Path,
    session: Session,
    run: db_orm.Run,
    fasta_dir: Path,
    hash_to_filename: dict[str, str],
    filename_to_hash: dict[str, str],
    query_hashes: list[str],
    subject_hash: str,
    *,
    quiet: bool = False,
) -> int:
    """Run fastANI many-vs-subject and log column of comparisons to database."""
    uname = platform.uname()
    uname_system = uname.system
    uname_release = uname.release
    uname_machine = uname.machine

    tool = tools.get_fastani()
    _check_tool_version(tool, run.configuration)

    config_id = run.configuration.configuration_id
    fragsize = run.configuration.fragsize
    if not fragsize:
        msg = f"ERROR: fastANI run-id {run.run_id} is missing fragsize parameter"
        sys.exit(msg)
    kmersize = run.configuration.kmersize
    if not kmersize:
        msg = f"ERROR: fastANI run-id {run.run_id} is missing kmersize parameter"
        sys.exit(msg)
    minmatch = run.configuration.minmatch
    if not minmatch:
        msg = f"ERROR: fastANI run-id {run.run_id} is missing minmatch parameter"
        sys.exit(msg)

    tmp_output = tmp_dir / f"queries_vs_{subject_hash}.csv"
    tmp_queries = tmp_dir / f"queries_vs_{subject_hash}.txt"
    with tmp_queries.open("w") as handle:
        for query_hash in query_hashes:
            handle.write(f"{fasta_dir / hash_to_filename[query_hash]}\n")

    if not quiet:
        print(
            f"INFO: Calling fastANI for {len(query_hashes)} queries vs {subject_hash}"
        )

    check_output(
        [
            str(tool.exe_path),
            "--ql",
            str(tmp_queries),
            "-r",
            str(fasta_dir / hash_to_filename[subject_hash]),
            # Send to file or just capture stdout?
            "-o",
            str(tmp_output),
            "--fragLen",
            str(fragsize),
            "-k",
            str(kmersize),
            "--minFraction",
            str(minmatch),
        ],
    )

    # Now do a bulk import. There shouldn't be any unless perhaps we have a
    # race condition with another thread, but skip any pre-existing entries via
    # Sqlite3's "INSERT OR IGNORE" using dialect's on_conflict_do_nothing method.
    session.execute(
        sqlite_insert(db_orm.Comparison).on_conflict_do_nothing(),
        [
            {
                "query_hash": query_hash,
                "subject_hash": subject_hash,
                "identity": identity,
                # Proxy values:
                "aln_length": None
                if orthologous_matches is None
                else round(run.configuration.fragsize * orthologous_matches),
                "sim_errors": None
                if fragments is None or orthologous_matches is None
                else fragments - orthologous_matches,
                "cov_query": None
                if fragments is None or orthologous_matches is None
                else orthologous_matches / fragments,
                "configuration_id": config_id,
                "uname_system": uname_system,
                "uname_release": uname_release,
                "uname_machine": uname_machine,
            }
            for (
                query_hash,
                subject_hash,
                identity,
                orthologous_matches,
                fragments,
            ) in fastani.parse_fastani_file(
                tmp_output,
                filename_to_hash,
                # This is used to infer failed alignments:
                expected_pairs={(_, subject_hash) for _ in query_hashes},
            )
        ],
    )
    session.commit()
    return 0


def compute_anim(  # noqa: C901, PLR0913, PLR0915
    tmp_dir: Path,
    session: Session,
    run: db_orm.Run,
    fasta_dir: Path,
    hash_to_filename: dict[str, str],
    filename_to_hash: dict[str, str],  # noqa: ARG001
    query_hashes: list[str],
    subject_hash: str,
    *,
    quiet: bool = False,
) -> int:
    """Run ANIm many-vs-subject and log column of comparisons to database."""
    uname = platform.uname()
    uname_system = uname.system
    uname_release = uname.release
    uname_machine = uname.machine

    nucmer = tools.get_nucmer()
    delta_filter = tools.get_delta_filter()
    _check_tool_version(nucmer, run.configuration)

    config_id = run.configuration.configuration_id
    mode = run.configuration.mode
    if not mode:
        msg = f"ERROR: ANIm run-id {run.run_id} is missing mode parameter"
        sys.exit(msg)

    subject_length = (
        session.query(db_orm.Genome)
        .where(db_orm.Genome.genome_hash == subject_hash)
        .one()
        .length
    )

    # nucmer does not handle spaces in filenames, neither quoted nor
    # escaped as slash-space. Therefore symlink or decompress to <MD5>.fasta:
    subject_fasta = tmp_dir / f"{subject_hash}.fasta"
    stage_file(fasta_dir / hash_to_filename[subject_hash], subject_fasta)

    db_entries = []
    try:
        for query_hash in query_hashes:
            query_length = (
                session.query(db_orm.Genome)
                .where(db_orm.Genome.genome_hash == query_hash)
                .one()
                .length
            )
            if query_hash != subject_hash:
                # Another thread may create/delete that FASTA name for our query
                # - so make a unique name for the temp file:
                query_fasta = tmp_dir / f"{query_hash}_vs_{subject_hash}.fasta"
                stage_file(fasta_dir / hash_to_filename[query_hash], query_fasta)
            else:
                # Can reuse the subject's decompressed file/symlink
                query_fasta = subject_fasta

            stem = tmp_dir / f"{query_hash}_vs_{subject_hash}"
            delta = tmp_dir / f"{query_hash}_vs_{subject_hash}.delta"
            deltafilter = tmp_dir / f"{query_hash}_vs_{subject_hash}.filter"

            if not quiet:
                print(
                    f"INFO: Calling nucmer for {hash_to_filename[query_hash]}"
                    f" vs {hash_to_filename[subject_hash]}"
                )

            # Here mode will be "mum" (default) or "maxmatch", meaning nucmer --mum etc.
            check_output(
                [
                    str(nucmer.exe_path),
                    "-p",
                    str(stem),
                    f"--{mode}",
                    # subject THEN query:
                    str(subject_fasta),
                    str(query_fasta),
                ],
            )
            if not delta.is_file():
                msg = f"ERROR: nucmer didn't make {delta}"  # pragma: no cover
                sys.exit(msg)  # pragma: no cover

            if not quiet:
                print(
                    f"INFO: Calling delta filter for {hash_to_filename[query_hash]}"
                    f" vs {hash_to_filename[subject_hash]}"
                )

            # The constant -1 option is used for 1-to-1 alignments in the delta-filter,
            # with no other options available for the end user.
            output = check_output(
                [
                    str(delta_filter.exe_path),
                    "-1",
                    str(delta),
                ],
            )
            # Don't really need to write this to disk except to help with testing intermediates
            with deltafilter.open("w") as handle:
                handle.write(output)

            query_aligned_bases, subject_aligned_bases, identity, sim_errors = (
                anim.parse_delta(deltafilter)
            )

            db_entries.append(
                {
                    "query_hash": query_hash,
                    "subject_hash": subject_hash,
                    "identity": identity,
                    "aln_length": query_aligned_bases,
                    "sim_errors": sim_errors,
                    "cov_query": None
                    if query_aligned_bases is None
                    else float(query_aligned_bases) / query_length,
                    "cov_subject": None
                    if subject_aligned_bases is None
                    else float(subject_aligned_bases) / subject_length,
                    "configuration_id": config_id,
                    "uname_system": uname_system,
                    "uname_release": uname_release,
                    "uname_machine": uname_machine,
                }
            )

            # Waiting to log the whole column does reduce DB contention and
            # seems to avoid locking problems, but also means zero feedback.
            # Logging every 25 entries or any similar fixed size seems likely
            # to risk locking. The following should result in staggered commits:
            if query_hash == subject_hash:
                db_entries = db_orm.attempt_insert(
                    session, db_entries, db_orm.Comparison
                )
            elif hash_to_filename[query_hash].endswith(".gz"):
                query_fasta.unlink()  # remove our decompressed copy

    except KeyboardInterrupt:
        # Try to abort gracefully without wasting the work done.
        msg = f"Interrupted, will attempt to log {len(db_entries)} completed comparisons\n"
        sys.stderr.write(msg)
        run.status = "Worker interrupted"

    if db_entries:
        session.execute(
            sqlite_insert(db_orm.Comparison).on_conflict_do_nothing(), db_entries
        )
        session.commit()

    if hash_to_filename[subject_hash].endswith(".gz"):
        subject_fasta.unlink()  # remove our decompressed copy

    return 0


def compute_anib(  # noqa: PLR0913
    tmp_dir: Path,
    session: Session,
    run: db_orm.Run,
    fasta_dir: Path,
    hash_to_filename: dict[str, str],
    filename_to_hash: dict[str, str],  # noqa: ARG001
    query_hashes: list[str],
    subject_hash: str,
    *,
    quiet: bool = False,
) -> int:
    """Run ANIb many-vs-subject and log column of comparisons to database."""
    uname = platform.uname()
    uname_system = uname.system
    uname_release = uname.release
    uname_machine = uname.machine

    tool = tools.get_blastn()
    _check_tool_version(tool, run.configuration)

    config_id = run.configuration_id
    fragsize = run.configuration.fragsize
    if not fragsize:
        msg = f"ERROR: ANIb run-id {run.run_id} is missing fragsize parameter"
        sys.exit(msg)
    subject_length = (
        session.query(db_orm.Genome)
        .where(db_orm.Genome.genome_hash == subject_hash)
        .one()
        .length
    )
    outfmt = "6 " + " ".join(anib.BLAST_COLUMNS)

    # makeblastdb does not handle spaces in filenames, neither quoted nor
    # escaped as slash-space. Therefore symlink or decompress to <MD5>.fasta:
    subject_fasta = tmp_dir / f"{subject_hash}.fasta"
    stage_file(fasta_dir / hash_to_filename[subject_hash], subject_fasta)

    tmp_db = tmp_dir / subject_hash  # prefix for BLAST DB

    if not quiet:
        print(f"INFO: Calling makeblastdb for {hash_to_filename[subject_hash]}")

    check_output(
        [
            str(tools.get_makeblastdb().exe_path),
            "-in",
            str(subject_fasta),
            "-input_type",
            "fasta",
            "-dbtype",
            "nucl",
            "-title",
            subject_hash,
            "-out",
            str(tmp_db),
        ],
    )

    if hash_to_filename[subject_hash].endswith(".gz"):
        subject_fasta.unlink()  # remove our decompressed copy

    db_entries = []
    try:
        for query_hash in query_hashes:
            tmp_tsv = tmp_dir / f"{query_hash}_vs_{subject_hash}.tsv"

            # Potential race condition if other columns are being computed with the
            # same tmp_dir - so give the fragments file a unique name using PID:
            tmp_frag_query = (
                tmp_dir / f"{query_hash}-fragments-{fragsize}-pid{os.getpid()}.fna"
            )

            anib.fragment_fasta_file(
                fasta_dir / hash_to_filename[query_hash],
                tmp_frag_query,
                fragsize,
            )

            if not quiet:
                print(
                    f"INFO: Calling blastn for {hash_to_filename[query_hash]}"
                    f" vs {hash_to_filename[subject_hash]}"
                )

            check_output(
                [
                    str(tool.exe_path),
                    "-query",
                    str(tmp_frag_query),
                    "-db",
                    str(tmp_db),
                    "-out",
                    str(tmp_tsv),
                    "-task",
                    "blastn",
                    "-outfmt",
                    outfmt,
                    "-xdrop_gap_final",
                    "150",
                    "-dust",
                    "no",
                    "-evalue",
                    "1e-15",
                ],
            )

            identity, aln_length, sim_errors = anib.parse_blastn_file(tmp_tsv)

            query_length = (
                session.query(db_orm.Genome)
                .where(db_orm.Genome.genome_hash == query_hash)
                .one()
                .length
            )

            db_entries.append(
                {
                    "query_hash": query_hash,
                    "subject_hash": subject_hash,
                    "identity": identity,
                    "aln_length": aln_length,
                    "sim_errors": sim_errors,
                    "cov_query": None
                    if aln_length is None
                    else aln_length / query_length,
                    "cov_subject": None
                    if aln_length is None
                    else aln_length / subject_length,
                    "configuration_id": config_id,
                    "uname_system": uname_system,
                    "uname_release": uname_release,
                    "uname_machine": uname_machine,
                }
            )

            # Waiting to log the whole column does reduce DB contention and
            # seems to avoid locking problems, but also means zero feedback.
            # Logging every 25 entries or any similar fixed size seems likely
            # to risk locking. The following should result in staggered commits:
            if query_hash == subject_hash:
                db_entries = db_orm.attempt_insert(
                    session, db_entries, db_orm.Comparison
                )

    except KeyboardInterrupt:
        # Try to abort gracefully without wasting the work done.
        msg = f"Interrupted, will attempt to log {len(db_entries)} completed comparisons\n"
        sys.stderr.write(msg)
        run.status = "Worker interrupted"

    if db_entries:
        session.execute(
            sqlite_insert(db_orm.Comparison).on_conflict_do_nothing(), db_entries
        )
        session.commit()

    return 0


def compute_dnadiff(  # noqa: C901, PLR0912, PLR0913, PLR0915
    tmp_dir: Path,
    session: Session,
    run: db_orm.Run,
    fasta_dir: Path,
    hash_to_filename: dict[str, str],
    filename_to_hash: dict[str, str],  # noqa: ARG001
    query_hashes: list[str],
    subject_hash: str,
    *,
    quiet: bool = False,
) -> int:
    """Run dnadiff many-vs-subject and log column of comparisons to database."""
    uname = platform.uname()
    uname_system = uname.system
    uname_release = uname.release
    uname_machine = uname.machine

    nucmer = tools.get_nucmer()
    delta_filter = tools.get_delta_filter()
    show_diff = tools.get_show_diff()
    show_coords = tools.get_show_coords()
    _check_tool_version(nucmer, run.configuration)

    config_id = run.configuration.configuration_id

    # nucmer does not handle spaces in filenames, neither quoted nor
    # escaped as slash-space. Therefore symlink or decompress to <MD5>.fasta:
    subject_fasta = tmp_dir / f"{subject_hash}.fasta"
    stage_file(fasta_dir / hash_to_filename[subject_hash], subject_fasta)

    db_entries = []
    try:
        for query_hash in query_hashes:
            query_length = (
                session.query(db_orm.Genome)
                .where(db_orm.Genome.genome_hash == query_hash)
                .one()
                .length
            )
            if query_hash != subject_hash:
                # Another thread may create/delete that FASTA name for our query
                # - so make a unique name for the temp file:
                query_fasta = tmp_dir / f"{query_hash}_vs_{subject_hash}.fasta"
                stage_file(fasta_dir / hash_to_filename[query_hash], query_fasta)
            else:
                # Can reuse the subject's decompressed file/symlink
                query_fasta = subject_fasta

            stem = tmp_dir / f"{query_hash}_vs_{subject_hash}"
            delta = tmp_dir / f"{query_hash}_vs_{subject_hash}.delta"
            deltafilter = tmp_dir / f"{query_hash}_vs_{subject_hash}.filter"
            qdiff = tmp_dir / f"{query_hash}_vs_{subject_hash}.qdiff"
            mcoords = tmp_dir / f"{query_hash}_vs_{subject_hash}.mcoords"

            if not quiet:
                print(
                    f"INFO: Calling nucmer for {hash_to_filename[query_hash]}"
                    f" vs {hash_to_filename[subject_hash]}"
                )
            # This should not be run in the same tmp_dir as ANIm, as the nucmer output will clash
            check_output(
                [
                    str(nucmer.exe_path),
                    "-p",
                    str(stem),
                    "--maxmatch",
                    # subject THEN query:
                    str(subject_fasta),
                    str(query_fasta),
                ],
            )
            if not delta.is_file():
                msg = f"ERROR: nucmer didn't make {delta}"  # pragma: no cover
                sys.exit(msg)  # pragma: no cover

            if not quiet:
                print(
                    f"INFO: Calling delta-filter for {hash_to_filename[query_hash]}"
                    f" vs {hash_to_filename[subject_hash]}"
                )
            output = check_output(
                [
                    str(delta_filter.exe_path),
                    "-m",
                    str(delta),
                ],
            )
            # May be able to avoid writing this to disk, but helps with testing intermediates
            with deltafilter.open("w") as handle:
                handle.write(output)

            if not quiet:
                print(
                    f"INFO: Calling show-diff for {hash_to_filename[query_hash]}"
                    f" vs {hash_to_filename[subject_hash]}"
                )
            output = check_output(
                [
                    str(show_diff.exe_path),
                    "-qH",
                    str(deltafilter),
                ],
            )
            # May be able to avoid writing this to disk, but helps with testing intermediates
            with qdiff.open("w") as handle:
                handle.write(output)

            if not quiet:
                print(
                    f"INFO: Calling show-coords for {hash_to_filename[query_hash]}"
                    f" vs {hash_to_filename[subject_hash]}"
                )
            output = check_output(
                [
                    str(show_coords.exe_path),
                    "-rclTH",
                    str(deltafilter),
                ],
            )
            # May be able to avoid writing this to disk, but helps with testing intermediates
            with mcoords.open("w") as handle:
                handle.write(output)

            identity, aligned_bases_with_gaps = dnadiff.parse_mcoords(mcoords)
            gap_lengths = dnadiff.parse_qdiff(qdiff)

            # For comparisons of closely related genomes, qdiff files might
            # be empty as there are no gaps in the alignments. In this case, we
            # want to treat gap_lengths as 0. In cases of comparisons
            # of distantly related genomes, we report gap_lengths as None.
            aln_length = (
                None
                if gap_lengths is None and aligned_bases_with_gaps is None
                else (aligned_bases_with_gaps or 0) - (gap_lengths or 0)
            )
            sim_errors = (
                None
                if identity is None or aligned_bases_with_gaps is None
                else round(
                    ((aligned_bases_with_gaps or 0) - (gap_lengths or 0))
                    * (1 - identity)
                )
            )
            cov_query = (
                None
                if aligned_bases_with_gaps is None or query_length == 0
                else ((aligned_bases_with_gaps or 0) - (gap_lengths or 0))
                / query_length
            )

            db_entries.append(
                {
                    "query_hash": query_hash,
                    "subject_hash": subject_hash,
                    "identity": identity,
                    "aln_length": aln_length,
                    "sim_errors": sim_errors,
                    "cov_query": cov_query,
                    "cov_subject": None,
                    "configuration_id": config_id,
                    "uname_system": uname_system,
                    "uname_release": uname_release,
                    "uname_machine": uname_machine,
                }
            )

            # Waiting to log the whole column does reduce DB contention and
            # seems to avoid locking problems, but also means zero feedback.
            # Logging every 25 entries or any similar fixed size seems likely
            # to risk locking. The following should result in staggered commits:
            if query_hash == subject_hash:
                db_entries = db_orm.attempt_insert(
                    session, db_entries, db_orm.Comparison
                )
            elif hash_to_filename[query_hash].endswith(".gz"):
                query_fasta.unlink()  # remove our decompressed copy
    except KeyboardInterrupt:
        # Try to abort gracefully without wasting the work done.
        msg = f"Interrupted, will attempt to log {len(db_entries)} completed comparisons\n"
        sys.stderr.write(msg)
        run.status = "Worker interrupted"

    if db_entries:
        session.execute(
            sqlite_insert(db_orm.Comparison).on_conflict_do_nothing(), db_entries
        )
        session.commit()

    if hash_to_filename[subject_hash].endswith(".gz"):
        subject_fasta.unlink()  # remove our decompressed copy

    return 0


@app.command(rich_help_panel="Method specific logging")
def log_sourmash(
    database: REQ_ARG_TYPE_DATABASE,
    run_id: REQ_ARG_TYPE_RUN_ID,
    manysearch: Annotated[
        Path,
        typer.Option(
            help="Sourmash-plugin-branchwater manysearch CSV output file",
            show_default=False,
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ],
    *,
    quiet: OPT_ARG_TYPE_QUIET = False,
) -> int:
    """Log an all-vs-all sourmash pairwise comparison to database."""
    if database != ":memory:" and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    uname = platform.uname()
    uname_system = uname.system
    uname_release = uname.release
    uname_machine = uname.machine

    if not quiet:
        print(f"Logging sourmash to {database}")
    session = db_orm.connect_to_db(database)
    run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()
    if run.configuration.method != "sourmash":
        msg = f"ERROR: Run-id {run_id} expected {run.configuration.method} results"
        sys.exit(msg)

    _check_tool_version(tools.get_sourmash(), run.configuration)

    config_id = run.configuration.configuration_id
    filename_to_hash = {_.fasta_filename: _.genome_hash for _ in run.fasta_hashes}

    # Now do a bulk import... but must skip any pre-existing entries
    # otherwise would hit sqlite3.IntegrityError for breaking uniqueness!
    # Do this via the Sqlite3 supported SQL command "INSERT OR IGNORE"
    # using the dialect's on_conflict_do_nothing method.
    # Repeating those calculations is a waste, but a performance trade off
    session.execute(
        sqlite_insert(db_orm.Comparison).on_conflict_do_nothing(),
        [
            {
                "query_hash": query_hash,
                "subject_hash": subject_hash,
                "identity": max_containment,
                "cov_query": query_containment,
                "configuration_id": config_id,
                "uname_system": uname_system,
                "uname_release": uname_release,
                "uname_machine": uname_machine,
            }
            for (
                query_hash,
                subject_hash,
                query_containment,
                max_containment,
            ) in sourmash.parse_sourmash_manysearch_csv(
                manysearch,
                filename_to_hash,
                # This is used to infer failed alignments:
                expected_pairs={
                    (q, s)
                    for q in filename_to_hash.values()
                    for s in filename_to_hash.values()
                },
            )
        ],
    )

    session.commit()
    return 0


def compute_external_alignment(  # noqa: C901, PLR0912, PLR0913, PLR0915
    tmp_dir: Path,  # noqa: ARG001
    session: Session,
    run: db_orm.Run,
    fasta_dir: Path,  # noqa: ARG001
    hash_to_filename: dict[str, str],  # noqa: ARG001
    filename_to_hash: dict[str, str],  # noqa: ARG001
    query_hashes: list[str],
    subject_hash: str,
    *,
    quiet: bool = False,
) -> int:
    """Compute and log column of comparisons from given MSA to database.

    Will only look at query in query_hashes vs subject_hash, but will also
    record reciprocal comparison as this method is symmetric.
    """
    uname = platform.uname()
    uname_system = uname.system
    uname_release = uname.release
    uname_machine = uname.machine

    config_id = run.configuration.configuration_id
    if run.configuration.method != "external-alignment":
        msg = f"ERROR: Run-id {run.run_id} expected {run.configuration.method} results"
        sys.exit(msg)
    if run.configuration.program:
        msg = f"ERROR: configuration.program={run.configuration.program!r} unexpected"
        sys.exit(msg)
    if run.configuration.version:
        msg = f"ERROR: configuration.version={run.configuration.version!r} unexpected"
        sys.exit(msg)
    if not run.configuration.extra:
        msg = "ERROR: Missing configuration.extra setting"
        sys.exit(msg)
    args = dict(_.split("=", 1) for _ in run.configuration.extra.split(";", 2))
    if list(args) != ["md5", "label", "alignment"]:
        msg = f"ERROR: configuration.extra={run.configuration.extra!r} unexpected"
        sys.exit(msg)

    alignment = Path(args["alignment"])
    if not alignment.is_absolute():
        # If not absolute, assume MSA path relative to the DB.
        # Get the DB filename via the session connection binding
        url = str(session.bind.url)
        if not url.startswith("sqlite:///"):
            msg = (  # pragma: nocover
                f"Expected SQLite3 URL to start sqlite:/// not {url}"
            )
            raise ValueError(msg)  # pragma: nocover
        if not quiet:
            print(f"DEBUG: Treating {alignment} as relative to {url}")
        alignment = Path(url[10:]).parent / alignment

    md5 = args["md5"]
    label = args["label"]
    del args

    if not quiet:
        print(f"INFO: Parsing {alignment} (MD5={md5}, label={label})")
    if not alignment.is_file():
        msg = f"ERROR: Missing alignment file {alignment}"
        sys.exit(msg)
    if md5 != file_md5sum(alignment):
        msg = f"ERROR: MD5 checksum of {alignment} didn't match."
        sys.exit(msg)

    if label == "md5":
        mapping = lambda x: x  # noqa: E731
    elif label == "filename":
        mapping = {_.fasta_filename: _.genome_hash for _ in run.fasta_hashes}.get
    else:
        mapping = {
            filename_stem(_.fasta_filename): _.genome_hash for _ in run.fasta_hashes
        }.get

    import numpy as np  # lazy import, although might be implicitly loaded already?
    from Bio.SeqIO.FastaIO import SimpleFastaParser  # deliberate lazy import

    # First col, computes N, logs 2N-1 (A-vs-A, B-vs-A, C-vs-A, ..., Z-vs-A and mirrors)
    # Second col, computes N-1, logs 2N-3 (skips A-vs-B, computes B-vs-B, ..., Z-vs-B)
    # ...
    # N-th col, computes 1, logs 1 (skips A-vs-Z, ..., Y-vs-Z, computes Z-vs-Z only)
    #
    # Hopefully snakemake will submit jobs in that order which would be most efficient.
    #
    # Potentially we could instead break this into N/2 jobs (even) or (N+1)/2 jobs (odd)
    # Thus when N is even, the first job would consisting of columns 1 and N, next
    # job columns 2 and N-1, etc. Each job does N computations and records 2N comparisons.
    # Similarly when N is odd, although there we have a final half-sized single-column job.

    # We could interpret the column number as the MSA ordering, but our internal
    # API by subject hash - and we can't assume the MSA is in any particular order.
    # Easiest way to solve this is two linear scans of the file, with a seek(0)
    with alignment.open() as handle:
        subject_seq = subject_title = ""  # placeholder values
        s_non_gaps = subject_seq_gaps = None  # placeholder values
        # Should load the file in binary mode as will be working in bytes
        for query_title, query_seq in SimpleFastaParser(handle):
            query_hash = mapping(query_title.split(None, 1)[0])
            if not query_hash:
                msg = f"ERROR: Could not map {query_title.split(None, 1)[0]} as {label}"
                sys.exit(msg)
            if query_hash == subject_hash:
                # for use in rest of the loop - an array of bytes!
                subject_seq = query_seq
                s_array = np.array(list(subject_seq), "S1")
                s_non_gaps = s_array != b"-"
                subject_seq_gaps = subject_seq.count("-")
                subject_title = query_title  # for use in logging
                break
        else:
            msg = f"Did not find subject {subject_hash} in {alignment.name}"
            raise ValueError(msg)

        db_entries = []

        handle.seek(0)
        for query_title, query_seq in SimpleFastaParser(handle):
            query_hash = mapping(query_title.split(None, 1)[0])
            if query_hash < subject_hash or query_hash not in query_hashes:
                # Exploiting symmetry to avoid double computation,
                # or not asked to compute this pairing (as already in the DB)
                continue
            if query_hash == subject_hash:
                # 100% identity and coverage, but need to calculate aln_length
                db_entries.append(
                    {
                        "query_hash": query_hash,
                        "subject_hash": subject_hash,
                        "identity": 1.0,
                        "aln_length": len(query_seq) - query_seq.count("-"),
                        "sim_errors": 0,
                        "cov_query": 1.0,
                        "cov_subject": 1.0,
                        "configuration_id": config_id,
                        "uname_system": uname_system,
                        "uname_release": uname_release,
                        "uname_machine": uname_machine,
                    }
                )
            else:
                # Full calculation required
                if len(query_seq) != len(subject_seq):
                    msg = (
                        "ERROR: Bad external-alignment, different lengths"
                        f" {len(query_seq)} and {len(subject_seq)}"
                        f" from {query_title.split(None, 1)[0]}"
                        f" and {subject_title.split(None, 1)[0]}\n"
                    )
                    sys.exit(msg)

                q_array = np.array(list(query_seq), "S1")
                q_non_gaps = q_array != b"-"
                # & is AND
                # | is OR
                # ^ is XOR
                # ~ is NOT
                # e.g. ~(q_gaps | s_gaps) would be entries with no gaps
                naive_matches = q_array == s_array  # includes double gaps!
                matches = int((naive_matches & q_non_gaps).sum())
                one_gapped = q_non_gaps ^ s_non_gaps
                non_gap_mismatches = int((~naive_matches & ~one_gapped).sum())
                either_gapped = int(one_gapped.sum())
                del naive_matches, q_non_gaps, q_array

                # Now compute the alignment metrics from that
                query_cov = (matches + non_gap_mismatches) / (
                    len(query_seq) - query_seq.count("-")
                )
                subject_cov = (matches + non_gap_mismatches) / (
                    len(subject_seq) - subject_seq_gaps
                )
                aln_length = matches + non_gap_mismatches + either_gapped
                sim_errors = non_gap_mismatches + either_gapped

                db_entries.append(
                    {
                        "query_hash": query_hash,
                        "subject_hash": subject_hash,
                        "identity": matches / aln_length,
                        "aln_length": aln_length,
                        "sim_errors": sim_errors,
                        "cov_query": query_cov,
                        "cov_subject": subject_cov,
                        "configuration_id": config_id,
                        "uname_system": uname_system,
                        "uname_release": uname_release,
                        "uname_machine": uname_machine,
                    }
                )
                # Fill in the symmetric entry
                db_entries.append(
                    {
                        "query_hash": subject_hash,
                        "subject_hash": query_hash,
                        "identity": matches / aln_length,
                        "aln_length": aln_length,
                        "sim_errors": sim_errors,
                        "cov_query": subject_cov,
                        "cov_subject": query_cov,
                        "configuration_id": config_id,
                        "uname_system": uname_system,
                        "uname_release": uname_release,
                        "uname_machine": uname_machine,
                    }
                )
    if not quiet:
        print(f"DEBUG: Logging {len(db_entries)} comparisons vs {subject_title}")
    # Now do a bulk import... but must skip any pre-existing entries
    # otherwise would hit sqlite3.IntegrityError for breaking uniqueness!
    # Do this via the Sqlite3 supported SQL command "INSERT OR IGNORE"
    # using the dialect's on_conflict_do_nothing method.
    # Repeating those calculations is a waste, could skip if partially done?
    session.execute(
        sqlite_insert(db_orm.Comparison).on_conflict_do_nothing(),
        db_entries,
    )
    session.commit()
    return 0


if __name__ == "__main__":
    sys.exit(app())  # pragma: no cover
