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

import logging
import os
import platform
import signal
import sys
import tempfile
from contextlib import nullcontext
from pathlib import Path
from typing import Annotated

import typer
from sqlalchemy import func, select
from sqlalchemy.exc import NoResultFound
from sqlalchemy.orm import Session

from pyani_plus import LOG_FILE, db_orm, log_sys_exit, setup_logger, tools
from pyani_plus.public_cli_args import (
    NO_PATH,
    OPT_ARG_TYPE_CACHE,
    OPT_ARG_TYPE_CREATE_DB,
    OPT_ARG_TYPE_DEBUG,
    OPT_ARG_TYPE_LOG,
    OPT_ARG_TYPE_TEMP,
    OPT_ARG_TYPE_VERSION,
    REQ_ARG_TYPE_DATABASE,
    REQ_ARG_TYPE_FASTA_DIR,
)

ASCII_GAP = ord("-")  # 45

app = typer.Typer(
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)


# Typer callback adding version to the base command
@app.callback()
def common(  # noqa: D103
    ctx: typer.Context,
    *,
    version: OPT_ARG_TYPE_VERSION = False,
) -> None:
    pass  # pragma: no cover


REQ_ARG_TYPE_RUN_ID = Annotated[
    int,
    typer.Option("--run-id", "-r", help="Database run ID", show_default=False),
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
REQ_ARG_TYPE_JSON_OUT = Annotated[
    Path,
    typer.Option(
        help="Output JSON file", metavar="FILE", dir_okay=False, file_okay=True
    ),
]
REQ_ARG_TYPE_FASTA_FILES = Annotated[
    list[Path],
    typer.Argument(
        help="Path(s) to FASTA file(s)",
        metavar="FILE",
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
        metavar="FILE",
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
        metavar="FILE",
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


RECORDING_FAILED = 2  # return code for successful calculation but failed to save to DB


def _check_tool_version(
    logger: logging.Logger,
    tool: tools.ExternalToolData,
    configuration: db_orm.Configuration,
) -> None:
    """Confirm the tool and version matches the given configuration.

    >>> from pathlib import Path
    >>> from pyani_plus import setup_logger
    >>> from pyani_plus.tools import ExternalToolData
    >>> from pyani_plus.db_orm import Configuration
    >>> logger = setup_logger(None)
    >>> tool = ExternalToolData(Path("/bin/guestimator"), "1.3")
    >>> config = Configuration(method="guessing", program="guestimator", version="1.2")
    >>> _check_tool_version(logger, tool, config)
    Traceback (most recent call last):
    ...
    SystemExit: Run configuration was guestimator 1.2 but we have guestimator 1.3

    This will typically be used when resuming a run, to confirm the tool
    binary name (the stem) and version detected on the system match the
    configuration recorded for the run in the database.
    """
    if (
        configuration.program != tool.exe_path.stem
        or configuration.version != tool.version
    ):
        msg = (
            "Run configuration was"
            f" {configuration.program} {configuration.version}"
            f" but we have {tool.exe_path.stem} {tool.version}"
        )
        log_sys_exit(logger, msg)


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
    logger = setup_logger(None, terminal_level=logging.INFO)
    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"Database '{database}' does not exist, but not using --create-db"
        sys.exit(msg)

    msg = f"Logging configuration to '{database}'"
    logger.info(msg)
    session = db_orm.connect_to_db(logger, database)
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
    msg = f"Configuration identifier {config.configuration_id}"
    logger.info(msg)
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
    logger = setup_logger(None, terminal_level=logging.INFO)
    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"Database '{database}' does not exist, but not using --create-db"
        sys.exit(msg)

    msg = f"Logging genome to '{database}'"
    logger.info(msg)
    session = db_orm.connect_to_db(logger, database)

    from rich.progress import Progress

    from pyani_plus import PROGRESS_BAR_COLUMNS
    from pyani_plus.utils import file_md5sum

    file_total = 0
    if fasta:
        with Progress(*PROGRESS_BAR_COLUMNS) as progress:
            for filename in progress.track(fasta, description="Processing..."):
                file_total += 1
                md5 = file_md5sum(filename)
                db_orm.db_genome(logger, session, filename, md5, create=True)
    session.commit()
    session.close()
    msg = f"Processed {file_total} FASTA files"
    logger.info(msg)
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
    logger = setup_logger(None, terminal_level=logging.INFO)
    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"Database '{database}' does not exist, but not using --create-db"
        sys.exit(msg)

    msg = f"Logging run to '{database}'"
    logger.info(msg)
    session = db_orm.connect_to_db(logger, database)

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

    from rich.progress import Progress

    from pyani_plus import PROGRESS_BAR_COLUMNS
    from pyani_plus.utils import check_fasta, file_md5sum

    fasta_to_hash = {}
    fasta_names = check_fasta(logger, fasta)
    if fasta_names:
        # Reuse existing genome entries and/or log new ones
        with Progress(*PROGRESS_BAR_COLUMNS) as progress:
            for filename in progress.track(fasta_names, description="Processing..."):
                md5 = file_md5sum(filename)
                fasta_to_hash[filename] = md5
                db_orm.db_genome(logger, session, filename, md5, create=True)

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
    msg = f"Run identifier {run_id}"
    logger.info(msg)
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
    logger = setup_logger(None, terminal_level=logging.INFO)
    if database != ":memory:" and not Path(database).is_file():
        msg = f"Database '{database}' does not exist"
        sys.exit(msg)

    msg = f"Logging comparison to '{database}'"
    logger.info(msg)
    session = db_orm.connect_to_db(logger, database)
    # Give a better error message that if adding comparison fails:
    if (
        not session.query(db_orm.Configuration)
        .where(db_orm.Configuration.configuration_id == config_id)
        .count()
    ):
        msg = f"{database} does not contain configuration_id={config_id}"
        log_sys_exit(logger, msg)

    from pyani_plus.utils import file_md5sum

    query_md5 = file_md5sum(query_fasta)
    db_orm.db_genome(logger, session, query_fasta, query_md5)

    subject_md5 = file_md5sum(subject_fasta)
    db_orm.db_genome(logger, session, subject_fasta, subject_md5)

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


def export_json_db_entries(
    logger: logging.Logger,
    json_filename: Path,
    configuration: db_orm.Configuration,
    db_entries: list[dict[str, str | float | int | None]],
) -> None:
    """Serialise DB entries for recording in JSON for later import.

    The list of entries must all have the same configuration identifier as the given
    configuration, and the same uname values as the current machine. This is assumed.
    """
    config_dict = {
        "method": configuration.method,
        "program": configuration.program,
        "version": configuration.version,
        "fragsize": configuration.fragsize,
        "mode": configuration.mode,
        "kmersize": configuration.kmersize,
        "minmatch": configuration.minmatch,
        "extra": configuration.extra,
    }

    uname = platform.uname()  # this should be cached, so fast
    unwanted_keys = {
        "configuration_id",
        "uname_system",
        "uname_release",
        "uname_machine",
    }

    import json  # lazy

    serialised = json.dumps(
        {
            "configuration": config_dict,
            "uname": {
                "system": uname.system,
                "release": uname.release,
                "machine": uname.machine,
            },
            "comparisons": [
                {k: v for (k, v) in _.items() if k not in unwanted_keys}
                for _ in db_entries
            ],
        }
    )

    with json_filename.open("w") as handle:
        handle.write(serialised)
    msg = f"Saved {len(db_entries)} comparisons to {json_filename}"
    logger.debug(msg)


def import_json_comparisons(  # noqa: PLR0915
    logger: logging.Logger, session: Session, json_filename: Path
) -> int:
    """Import a JSON file of comparisons into the database."""
    import json  # lazy import

    msg = f"Importing {json_filename}"
    logger.debug(msg)

    # Consider incremental loading in case of large files?
    with json_filename.open("rb") as handle:
        raw = handle.read()
    if not raw:
        msg = f"JSON file '{json_filename}' is empty"
        logger.debug(msg)
        return 0
    try:
        data = json.loads(raw)
    except ValueError:
        logger.exception("Unable to parse JSON:")
        msg = f"JSON file '{json_filename}' invalid"
        log_sys_exit(logger, msg)
    del raw

    if (
        not isinstance(data, dict)
        or "configuration" not in data
        or not isinstance(data["configuration"], dict)
        or "uname" not in data
        or not isinstance(data["uname"], dict)
        or "comparisons" not in data
        or not isinstance(data["comparisons"], list)
    ):
        msg = f"JSON file '{json_filename}' does not use the expected structure"
        log_sys_exit(logger, msg)

    configuration = data["configuration"]
    comparisons = data["comparisons"]
    try:
        uname_system = data["uname"]["system"]
        uname_release = data["uname"]["release"]
        uname_machine = data["uname"]["machine"]
    except KeyError:
        msg = f"JSON file '{json_filename}' uname incomplete"
        session.close()
        log_sys_exit(logger, msg)
    del data

    try:
        config_id = db_orm.db_configuration(
            session=session,
            method=configuration["method"],
            program=configuration["program"],
            version=configuration["version"],
            fragsize=configuration.get("fragsize", None),
            mode=configuration.get("mode", None),
            kmersize=configuration.get("kmersize", None),
            minmatch=configuration.get("minmatch", None),
            extra=configuration.get("extra", None),
            create=False,
        ).configuration_id
    except KeyError:
        msg = f"JSON file '{json_filename}' configuration incomplete"
        session.close()
        log_sys_exit(logger, msg)
    except NoResultFound:
        msg = f"JSON file '{json_filename}' configuration not in database"
        session.close()
        log_sys_exit(logger, msg)
    del configuration
    msg = f"Configuration identifier {config_id} in database"
    logger.debug(msg)

    if not comparisons:
        msg = f"JSON file '{json_filename}' has no comparisons"
        session.close()
        logger.warning(msg)
        return 0

    try:
        db_entries = [
            {
                "query_hash": row["query_hash"],
                "subject_hash": row["subject_hash"],
                "identity": row["identity"],
                "aln_length": row.get("aln_length", None),
                "sim_errors": row.get("sim_errors", None),
                "cov_query": row.get("cov_query", None),
                "configuration_id": config_id,
                "uname_system": uname_system,
                "uname_release": uname_release,
                "uname_machine": uname_machine,
            }
            for row in comparisons
        ]
    except KeyError:
        msg = f"JSON file '{json_filename}' comparison(s) incomplete"
        session.close()
        log_sys_exit(logger, msg)

    if not db_orm.insert_comparisons_with_retries(
        logger, session, db_entries, source=str(json_filename)
    ):  # pragma: no cover
        msg = f"Failed to record '{json_filename}' comparisons to database"
        session.close()
        log_sys_exit(logger, msg)

    return len(comparisons)


@app.command(rich_help_panel="Low-level logging")
def import_comparisons(
    database: REQ_ARG_TYPE_DATABASE,
    json: Annotated[
        list[Path],
        typer.Argument(
            help="Path(s) of JSON file(s) containing comparisons.",
            metavar="FILE",
            show_default=False,
            exists=True,
            dir_okay=False,
            file_okay=True,
        ),
    ],
    *,
    debug: OPT_ARG_TYPE_DEBUG = False,
    log: OPT_ARG_TYPE_LOG = NO_PATH,
) -> int:
    """Parse one or more JSON comparison files and log them to the database.

    The database must already contain a matching configuration table entry
    (see the log-configuration command to do this at the command line) and
    matching genome entries (see the log-genome command).

    The JSON files will contain values for a configuration table entry, uname
    values applied to all the comparisons, and multiple entries to be recorded
    in the comparisons table.
    """
    logger = setup_logger(log, terminal_level=logging.DEBUG if debug else logging.ERROR)
    if database != ":memory:" and not Path(database).is_file():
        msg = f"Database '{database}' does not exist"
        sys.exit(msg)

    msg = f"Logging comparison to '{database}'"
    logger.info(msg)
    session = db_orm.connect_to_db(logger, database)

    if (
        session.execute(select(func.count()).select_from(db_orm.Configuration)).scalar()
        == 0
    ):
        msg = f"{database} does not contain any configurations"
        log_sys_exit(logger, msg)

    if session.execute(select(func.count()).select_from(db_orm.Genome)).scalar() == 0:
        msg = f"{database} does not contain any genomes"
        log_sys_exit(logger, msg)

    for filename in json:
        count = import_json_comparisons(logger, session, filename)
        msg = f"Imported {count} from '{filename}'"
        logger.info(msg)
    return 0


@app.command()
def prepare_genomes(
    database: REQ_ARG_TYPE_DATABASE,
    run_id: Annotated[
        int | None,
        typer.Option(help="Which run to prepare", show_default=False),
    ],
    cache: OPT_ARG_TYPE_CACHE = Path(),
    *,
    debug: OPT_ARG_TYPE_DEBUG = False,
    log: OPT_ARG_TYPE_LOG = NO_PATH,
) -> int:
    """Prepare any intermediate files needed prior to computing ANI values.

    This requires you already have a run defined in the database, and would
    be followed by running the private CLI ``compute-column`` command. Most
    methods do not need a prepare step, but for example ``sourmash`` does.
    """
    # Should this be splittable for running on the cluster? I assume most
    # cases this is IO bound rather than CPU bound so is this helpful?
    logger = setup_logger(
        log, terminal_level=logging.DEBUG if debug else logging.ERROR, plain=True
    )
    if database != ":memory:" and not Path(database).is_file():
        msg = f"Database '{database}' does not exist"
        log_sys_exit(logger, msg)
    session = db_orm.connect_to_db(logger, database)
    run = db_orm.load_run(session, run_id)
    try:
        return prepare(logger, run, cache)
    except Exception as err:  # pragma: nocover
        logger.exception("Unhandled exception:")
        msg = f"Unhandled exception preparing genomes: {err}"
        log_sys_exit(logger, msg)
        return 1  # for mypy


def prepare(logger: logging.Logger, run: db_orm.Run, cache: Path) -> int:
    """Call prepare-genomes with a progress bar."""
    n = run.genomes.count()
    done = run.comparisons().count()
    if done == n**2:
        msg = f"Skipping preparation, run already has all {n**2}={n}² pairwise values"
        logger.info(msg)
        return 0
    config = run.configuration
    method = config.method

    import importlib

    try:
        module = importlib.import_module(
            f"pyani_plus.methods.{method.lower().replace('-', '_')}"
        )
    except ModuleNotFoundError:
        msg = f"Unknown method {method}, check tool version?"
        log_sys_exit(logger, msg)
    if not hasattr(module, "prepare_genomes"):
        msg = f"No per-genome preparation required for {method}"
        logger.info(msg)  # debug level?
        return 0

    msg = f"Preparing {n} genomes under cache '{cache}'"
    logger.info(msg)

    from rich.progress import Progress

    from pyani_plus import PROGRESS_BAR_COLUMNS

    with Progress(*PROGRESS_BAR_COLUMNS) as progress:
        for _ in progress.track(
            module.prepare_genomes(logger, run, cache),
            description="Processing...  ",  # spaces to match "Indexing FASTAs" etc
            total=n,
        ):
            pass
    logger.debug("Done")
    return 0


@app.command()
def compute_column(  # noqa: C901, PLR0913, PLR0912, PLR0915
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
    json: REQ_ARG_TYPE_JSON_OUT,
    *,
    cache: OPT_ARG_TYPE_CACHE = Path(),
    temp: OPT_ARG_TYPE_TEMP = Path("-"),
    debug: OPT_ARG_TYPE_DEBUG = False,
    log: OPT_ARG_TYPE_LOG = LOG_FILE,
) -> int:
    """Run the method for one column and log pairwise comparisons to the database.

    The database and run ID specify the method and configuration, additionally
    you must supply a subject filename, hash, or column index to control which
    column of the matrix is to be computed.

    If using a column number, these are taken to be one based meaning in the range
    1 up to and including the number of genomes in the run. Using 0 instead means
    compute all the columns.

    Some methods like sourmash require you first run the prepare-genomes command
    (which for sourmash builds signature files from each FASTA file). You must
    use the same cache location for that and when you run compute-column.
    """
    # Do NOT write to the main thread's log (risk of race conditions appending
    # to the same file, locking, etc) - will use a column-specific log soon!
    try:
        column = int(subject)
    except ValueError:
        # Will have to wait and setup the file-based log later..
        column = -1
    if column < 0:
        log_ready = False
    else:
        log_ready = True
        log = Path(str(log)[: -len(log.suffix)] + f".{column}" + log.suffix)
    logger = setup_logger(
        log if log_ready else None,
        terminal_level=logging.DEBUG if debug else logging.ERROR,
        plain=True,
    )
    msg = f"Starting compute-column for {subject} to {json}"
    logger.debug(msg)

    if database != ":memory:" and not Path(database).is_file():
        msg = f"Database '{database}' does not exist"
        log_sys_exit(logger, msg)

    # We want to receive any SIGINT as a KeyboardInterrupt even if we
    # are run via a bash shell or other non-interactive setting:
    signal.signal(signal.SIGINT, signal.default_int_handler)

    # Under SLURM when a job is cancelled via scancel, it sends SIGTERM
    # for a graceful shutdown then 30s later a SIGKILL hard kill.
    # For simplicity, treat SIGTERM as a KeyboardInterrupt too.
    signal.signal(signal.SIGTERM, signal.default_int_handler)

    try:
        session = db_orm.connect_to_db(logger, database)
        run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()
        config = run.configuration
        method = config.method

        filename_to_hash = {_.fasta_filename: _.genome_hash for _ in run.fasta_hashes}
        hash_to_filename = {_.genome_hash: _.fasta_filename for _ in run.fasta_hashes}
        n = len(hash_to_filename)

        if subject in hash_to_filename:
            subject_hash = subject
            column = sorted(hash_to_filename).index(subject_hash) + 1
        elif Path(subject).name in filename_to_hash:
            subject_hash = filename_to_hash[Path(subject).name]
            column = sorted(hash_to_filename).index(subject_hash) + 1
        else:
            try:
                column = int(subject)
            except ValueError:
                msg = f"Did not recognise {subject!r} as an MD5 hash, filename, or column number in run-id {run_id}"
                log_sys_exit(logger, msg)
            if 0 < column <= n:
                subject_hash = sorted(hash_to_filename)[column - 1]
            elif column == 0:
                if method == "sourmash":
                    subject_hash = ""
                else:
                    msg = "All columns currently only implemented for sourmash"
                    log_sys_exit(logger, msg)
            else:
                msg = (
                    f"Single column should be in range 1 to {n},"
                    f" or for some methods {0} meaning all columns, but not {subject}"
                )
                log_sys_exit(logger, msg)

        if not log_ready:
            # Column worker specific log files!
            log = Path(str(log)[: -len(log.suffix)] + f".{column}" + log.suffix)
            logger = setup_logger(
                log,
                terminal_level=logging.DEBUG if debug else logging.ERROR,
                plain=True,
            )
        msg = f"Logging {method} compute-column {column} to {log}"
        logger.info(msg)

        if column == 0:
            # Computing all the matrix, but are there might be some rows/cols already done?
            query_hashes = {
                _.genome_hash: _.length for _ in run.genomes
            }  # assume all needed
        else:
            # What comparisons are needed? Record the query genome lengths too
            # (doing this once at the start to avoid a small lookup for each query)
            missing_query_hashes = set(hash_to_filename).difference(
                comp.query_hash
                for comp in run.comparisons().where(
                    db_orm.Comparison.subject_hash == subject_hash
                )
            )
            query_hashes = {
                _.genome_hash: _.length
                for _ in run.genomes
                if _.genome_hash in missing_query_hashes
            }
            del missing_query_hashes

        if not query_hashes:
            msg = f"No {method} comparisons needed against {subject_hash}"
            logger.info(msg)
            return 0

        # Will probably want to move each of these functions to the relevant method module...
        try:
            compute = {
                "fastANI": compute_fastani,
                "ANIb": compute_anib,
                "ANIm": compute_anim,
                "dnadiff": compute_dnadiff,
                "sourmash": compute_sourmash,
                "external-alignment": compute_external_alignment,
            }[method]
        except KeyError:
            msg = f"Unknown method {method} for run-id {run_id} in {database}"
            log_sys_exit(logger, msg)

        # On a cluster we are likely in a temp working directory, meaning
        # if it is a relative path, run.fasta_directory is useless without
        # evaluating it relative to the DB filename.
        fasta_dir = Path(run.fasta_directory)
        if not fasta_dir.is_absolute():
            fasta_dir = (database.parent / fasta_dir).absolute()
        msg = f"FASTA folder {fasta_dir}"
        logger.debug(msg)

        # Either use the specified temp-directory (and do not clean up),
        # or use a system temp-directory (and do clean up)
        if temp == Path("-"):  # dummy value accepted at command line
            temp = None
        if temp:
            temp = temp / f"c{column}"  # avoid worries about name clashes
            temp.mkdir(exist_ok=True)
            msg = f"Using temp folder {temp}"
            logger.debug(msg)

        msg = (
            f"Calling {method} for {len(query_hashes)} queries"
            if column == 0
            else f"Calling {method} for {len(query_hashes)} queries vs {subject_hash}."
        )
        logger.info(msg)
    except Exception as err:  # pragma: nocover
        logger.exception("Unhandled exception:")
        msg = f"Unhandled exception before compute: {err}"
        log_sys_exit(logger, msg)
        return 1  # for mypy

    try:
        with nullcontext(temp) if temp else tempfile.TemporaryDirectory() as tmp_dir:
            return compute(
                logger,
                Path(tmp_dir),
                session,
                run,
                json,
                fasta_dir,
                hash_to_filename,
                filename_to_hash,
                query_hashes,
                subject_hash,
                cache=cache,
            )
    except Exception as err:  # pragma: nocover
        logger.exception("Unhandled exception:")
        msg = f"Unhandled exception during compute: {err}"
        log_sys_exit(logger, msg)
        return 1  # for mypy


def compute_fastani(  # noqa: PLR0913, PLR0915
    logger: logging.Logger,
    tmp_dir: Path,
    session: Session,
    run: db_orm.Run,
    json_filename: Path,
    fasta_dir: Path,
    hash_to_filename: dict[str, str],
    filename_to_hash: dict[str, str],
    query_hashes: dict[str, int],
    subject_hash: str,
    *,
    cache: Path = Path(),  # noqa: ARG001
) -> int:
    """Run fastANI many-vs-subject and log column of comparisons to JSON.

    Unfortunately fastANI does not produce any output to file until it finishes
    all the given comparisons. This is problematic in terms of no user feedback,
    but also means if the job is interrupted there is nothing we can salvage.
    It also meant we tried to write all the comparisons in one batch at the end
    which means an all-or-nothing success/failure mode. Therefore run in batches.
    """
    uname = platform.uname()
    uname_system = uname.system
    uname_release = uname.release
    uname_machine = uname.machine

    configuration = run.configuration

    tool = tools.get_fastani()
    _check_tool_version(logger, tool, configuration)

    config_id = configuration.configuration_id
    fragsize = configuration.fragsize
    if not fragsize:
        msg = f"fastANI run-id {run.run_id} is missing fragsize parameter"
        log_sys_exit(logger, msg)
    kmersize = configuration.kmersize
    if not kmersize:
        msg = f"fastANI run-id {run.run_id} is missing kmersize parameter"
        log_sys_exit(logger, msg)
    minmatch = configuration.minmatch
    if not minmatch:
        msg = f"fastANI run-id {run.run_id} is missing minmatch parameter"
        log_sys_exit(logger, msg)

    from pyani_plus.methods import fastani  # lazy import
    from pyani_plus.utils import check_output

    # Given query_hashes as a dict (hash to query length), but only need hashes here:
    pending = sorted(query_hashes)
    del query_hashes
    db_entries: list[dict[str, str | float | None]] = []
    batch_size = 500
    try:
        while pending:
            batch = pending[:batch_size]
            pending = pending[batch_size:]

            # Could use a new filename for each batch?
            tmp_output = tmp_dir / f"queries_vs_{subject_hash}.csv"
            tmp_queries = tmp_dir / f"queries_vs_{subject_hash}.txt"
            with tmp_queries.open("w") as handle:
                for query_hash in batch:
                    handle.write(f"{fasta_dir / hash_to_filename[query_hash]}\n")

            msg = f"About to call fastANI on batch of {len(batch)} vs {subject_hash}"
            logger.debug(msg)
            check_output(
                logger,
                [
                    str(tool.exe_path),
                    "--ql",
                    str(tmp_queries),
                    "-r",
                    str(fasta_dir / hash_to_filename[subject_hash]),
                    # Send to file or just capture stdout?
                    # No point using stdout as only written after all computed
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
            msg = f"Called fastANI on batch of {len(batch)} vs {subject_hash}, about to parse output."
            logger.debug(msg)
            db_entries.extend(
                {
                    "query_hash": query_hash,
                    "subject_hash": subject_hash,
                    "identity": identity,
                    # Proxy values:
                    "aln_length": None
                    if orthologous_matches is None
                    else round(configuration.fragsize * orthologous_matches),
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
                    expected_pairs={(_, subject_hash) for _ in batch},
                )
            )
            msg = f"Parsed fastANI batch of {len(batch)} vs {subject_hash}; recording {len(db_entries)} to JSON."
            logger.debug(msg)
            export_json_db_entries(logger, json_filename, configuration, db_entries)

    except KeyboardInterrupt:  # pragma: no cover
        # Try to abort gracefully without wasting the work done.
        msg = f"Interrupted with {len(db_entries)} completed fastANI comparisons"
        logger.error(msg)  # noqa: TRY400
        run.status = "Worker interrupted"
        session.commit()

    try:
        # This is redundant in the good path...
        export_json_db_entries(logger, json_filename, configuration, db_entries)
    except Exception:  # pragma: no cover
        logger.exception("Unexpected exception saving JSON:")
        return RECORDING_FAILED
    else:
        return 0


def compute_anim(  # noqa: C901, PLR0913, PLR0915
    logger: logging.Logger,
    tmp_dir: Path,
    session: Session,
    run: db_orm.Run,
    json_filename: Path,
    fasta_dir: Path,
    hash_to_filename: dict[str, str],
    filename_to_hash: dict[str, str],  # noqa: ARG001
    query_hashes: dict[str, int],
    subject_hash: str,
    *,
    cache: Path = Path(),  # noqa: ARG001
) -> int:
    """Run ANIm many-vs-subject and log column of comparisons to JSON."""
    uname = platform.uname()
    uname_system = uname.system
    uname_release = uname.release
    uname_machine = uname.machine

    batch_size = 50

    configuration = run.configuration

    nucmer = tools.get_nucmer()
    delta_filter = tools.get_delta_filter()
    _check_tool_version(logger, nucmer, configuration)

    config_id = configuration.configuration_id
    mode = configuration.mode
    if not mode:
        msg = f"ANIm run-id {run.run_id} is missing mode parameter"
        log_sys_exit(logger, msg)

    subject_length = (
        session.query(db_orm.Genome)
        .where(db_orm.Genome.genome_hash == subject_hash)
        .one()
        .length
    )

    from pyani_plus.methods import anim  # lazy import
    from pyani_plus.utils import check_output, stage_file

    # nucmer does not handle spaces in filenames, neither quoted nor
    # escaped as slash-space. Therefore symlink or decompress to <MD5>.fasta:
    subject_fasta = tmp_dir / f"{subject_hash}.fasta"
    stage_file(logger, fasta_dir / hash_to_filename[subject_hash], subject_fasta)

    db_entries = []
    new = 0
    try:
        for query_hash, query_length in query_hashes.items():
            if query_hash != subject_hash:
                # Another thread may create/delete that FASTA name for our query
                # - so make a unique name for the temp file:
                query_fasta = tmp_dir / f"{query_hash}_vs_{subject_hash}.fasta"
                stage_file(
                    logger, fasta_dir / hash_to_filename[query_hash], query_fasta
                )
            else:
                # Can reuse the subject's decompressed file/symlink
                query_fasta = subject_fasta

            stem = tmp_dir / f"{query_hash}_vs_{subject_hash}"
            delta = tmp_dir / f"{query_hash}_vs_{subject_hash}.delta"
            deltafilter = tmp_dir / f"{query_hash}_vs_{subject_hash}.filter"

            msg = (
                f"Calling nucmer for {hash_to_filename[query_hash]}"
                f" vs {hash_to_filename[subject_hash]}"
            )
            logger.info(msg)

            # Here mode will be "mum" (default) or "maxmatch", meaning nucmer --mum etc.
            check_output(
                logger,
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
                msg = f"nucmer didn't make {delta}"  # pragma: no cover
                log_sys_exit(logger, msg)  # pragma: no cover

            msg = (
                f"Calling delta filter for {hash_to_filename[query_hash]}"
                f" vs {hash_to_filename[subject_hash]}"
            )
            logger.info(msg)
            # The constant -1 option is used for 1-to-1 alignments in the delta-filter,
            # with no other options available for the end user.
            output = check_output(
                logger,
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

            if query_hash != subject_hash and hash_to_filename[query_hash].endswith(
                ".gz"
            ):
                query_fasta.unlink()  # remove our decompressed copy

            new += 1
            if new >= batch_size:  # pragma: no cover
                export_json_db_entries(logger, json_filename, configuration, db_entries)
                new = 0
    except KeyboardInterrupt:
        # Try to abort gracefhave logged JSON file with  the work done.
        msg = f"Interrupted with {len(db_entries)} completed ANIm comparisons"
        logger.error(msg)  # noqa: TRY400
        run.status = "Worker interrupted"
        session.commit()

    if hash_to_filename[subject_hash].endswith(".gz"):
        subject_fasta.unlink()  # remove our decompressed copy

    if not new:  # pragma: no cover
        return 0
    try:
        export_json_db_entries(logger, json_filename, configuration, db_entries)
    except Exception:  # pragma: no cover
        logger.exception("Unexpected exception saving JSON:")
        return RECORDING_FAILED
    else:
        return 0


def compute_anib(  # noqa: PLR0913, PLR0915
    logger: logging.Logger,
    tmp_dir: Path,
    session: Session,
    run: db_orm.Run,
    json_filename: Path,
    fasta_dir: Path,
    hash_to_filename: dict[str, str],
    filename_to_hash: dict[str, str],  # noqa: ARG001
    query_hashes: dict[str, int],
    subject_hash: str,
    *,
    cache: Path = Path(),  # noqa: ARG001
) -> int:
    """Run ANIb many-vs-subject and log column of comparisons to JSON."""
    uname = platform.uname()
    uname_system = uname.system
    uname_release = uname.release
    uname_machine = uname.machine

    batch_size = 50

    configuration = run.configuration

    tool = tools.get_blastn()
    _check_tool_version(logger, tool, configuration)

    config_id = run.configuration_id
    fragsize = configuration.fragsize
    if not fragsize:
        msg = f"ANIb run-id {run.run_id} is missing fragsize parameter"
        log_sys_exit(logger, msg)
    subject_length = (
        session.query(db_orm.Genome)
        .where(db_orm.Genome.genome_hash == subject_hash)
        .one()
        .length
    )

    from pyani_plus.methods import anib  # lazy import
    from pyani_plus.utils import check_output, stage_file

    outfmt = "6 " + " ".join(anib.BLAST_COLUMNS)

    # makeblastdb does not handle spaces in filenames, neither quoted nor
    # escaped as slash-space. Therefore symlink or decompress to <MD5>.fasta:
    subject_fasta = tmp_dir / f"{subject_hash}.fasta"
    stage_file(logger, fasta_dir / hash_to_filename[subject_hash], subject_fasta)

    tmp_db = tmp_dir / subject_hash  # prefix for BLAST DB

    msg = f"Calling makeblastdb for {hash_to_filename[subject_hash]}"
    logger.info(msg)
    check_output(
        logger,
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
    new = 0
    try:
        for query_hash, query_length in query_hashes.items():
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

            msg = (
                f"Calling blastn for {hash_to_filename[query_hash]}"
                f" vs {hash_to_filename[subject_hash]}"
            )
            logger.info(msg)
            check_output(
                logger,
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
            new += 1
            if new >= batch_size:  # pragma: no cover
                export_json_db_entries(logger, json_filename, configuration, db_entries)
                new = 0

    except KeyboardInterrupt:
        # Try to abort gracefully without wasting the work done.
        msg = f"Interrupted with {len(db_entries)} completed ANIb comparisons"
        logger.error(msg)  # noqa: TRY400
        run.status = "Worker interrupted"
        session.commit()

    if not new:  # pragma: no cover
        return 0
    try:
        export_json_db_entries(logger, json_filename, configuration, db_entries)
    except Exception:  # pragma: no cover
        logger.exception("Unexpected exception saving JSON:")
        return RECORDING_FAILED
    else:
        return 0


def compute_dnadiff(  # noqa: PLR0913, PLR0915
    logger: logging.Logger,
    tmp_dir: Path,
    session: Session,
    run: db_orm.Run,
    json_filename: Path,
    fasta_dir: Path,
    hash_to_filename: dict[str, str],
    filename_to_hash: dict[str, str],  # noqa: ARG001
    query_hashes: dict[str, int],
    subject_hash: str,
    *,
    cache: Path = Path(),  # noqa: ARG001
) -> int:
    """Run dnadiff many-vs-subject and log column of comparisons to JSON."""
    uname = platform.uname()
    uname_system = uname.system
    uname_release = uname.release
    uname_machine = uname.machine

    batch_size = 50

    configuration = run.configuration

    nucmer = tools.get_nucmer()
    delta_filter = tools.get_delta_filter()
    show_diff = tools.get_show_diff()
    show_coords = tools.get_show_coords()
    _check_tool_version(logger, nucmer, configuration)

    config_id = configuration.configuration_id

    from pyani_plus.methods import dnadiff  # lazy import
    from pyani_plus.utils import check_output, stage_file

    # nucmer does not handle spaces in filenames, neither quoted nor
    # escaped as slash-space. Therefore symlink or decompress to <MD5>.fasta:
    subject_fasta = tmp_dir / f"{subject_hash}.fasta"
    stage_file(logger, fasta_dir / hash_to_filename[subject_hash], subject_fasta)

    db_entries = []
    new = 0
    try:
        for query_hash, query_length in query_hashes.items():
            if query_hash != subject_hash:
                # Another thread may create/delete that FASTA name for our query
                # - so make a unique name for the temp file:
                query_fasta = tmp_dir / f"{query_hash}_vs_{subject_hash}.fasta"
                stage_file(
                    logger, fasta_dir / hash_to_filename[query_hash], query_fasta
                )
            else:
                # Can reuse the subject's decompressed file/symlink
                query_fasta = subject_fasta

            stem = tmp_dir / f"{query_hash}_vs_{subject_hash}"
            delta = tmp_dir / f"{query_hash}_vs_{subject_hash}.delta"
            deltafilter = tmp_dir / f"{query_hash}_vs_{subject_hash}.filter"
            qdiff = tmp_dir / f"{query_hash}_vs_{subject_hash}.qdiff"
            mcoords = tmp_dir / f"{query_hash}_vs_{subject_hash}.mcoords"

            msg = (
                f"Calling nucmer for {hash_to_filename[query_hash]}"
                f" vs {hash_to_filename[subject_hash]}"
            )
            logger.info(msg)
            # This should not be run in the same tmp_dir as ANIm, as the nucmer output will clash
            check_output(
                logger,
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
                msg = f"nucmer didn't make {delta}"  # pragma: no cover
                log_sys_exit(logger, msg)  # pragma: no cover

            msg = (
                f"Calling delta-filter for {hash_to_filename[query_hash]}"
                f" vs {hash_to_filename[subject_hash]}"
            )
            logger.info(msg)
            output = check_output(
                logger,
                [
                    str(delta_filter.exe_path),
                    "-m",
                    str(delta),
                ],
            )
            # May be able to avoid writing this to disk, but helps with testing intermediates
            with deltafilter.open("w") as handle:
                handle.write(output)

            msg = (
                f"Calling show-diff for {hash_to_filename[query_hash]}"
                f" vs {hash_to_filename[subject_hash]}"
            )
            logger.info(msg)
            output = check_output(
                logger,
                [
                    str(show_diff.exe_path),
                    "-qH",
                    str(deltafilter),
                ],
            )
            # May be able to avoid writing this to disk, but helps with testing intermediates
            with qdiff.open("w") as handle:
                handle.write(output)

            msg = (
                f"Calling show-coords for {hash_to_filename[query_hash]}"
                f" vs {hash_to_filename[subject_hash]}"
            )
            logger.info(msg)
            output = check_output(
                logger,
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

            if query_hash != subject_hash and hash_to_filename[query_hash].endswith(
                ".gz"
            ):
                query_fasta.unlink()  # remove our decompressed copy

            new += 1
            if new >= batch_size:  # pragma: no cover
                export_json_db_entries(logger, json_filename, configuration, db_entries)
                new = 0

    except KeyboardInterrupt:
        # Try to abort gracefully without wasting the work done.
        msg = f"Interrupted with {len(db_entries)} completed dnadiff comparisons"
        logger.error(msg)  # noqa: TRY400
        run.status = "Worker interrupted"
        session.commit()

    if not new:  # pragma: no cover
        return 0
    try:
        export_json_db_entries(logger, json_filename, configuration, db_entries)
    except Exception:  # pragma: no cover
        logger.exception("Unexpected exception saving JSON:")
        return RECORDING_FAILED
    else:
        return 0


def compute_sourmash(  # noqa: PLR0913
    logger: logging.Logger,
    tmp_dir: Path,
    session: Session,
    run: db_orm.Run,
    json_filename: Path,
    fasta_dir: Path,  # noqa: ARG001
    hash_to_filename: dict[str, str],  # noqa: ARG001
    filename_to_hash: dict[str, str],  # noqa: ARG001
    query_hashes: dict[str, int],
    subject_hash: str,
    *,
    cache: Path = Path(),
) -> int:
    """Run many-vs-subject for sourmash and log to JSON."""
    # Not using try/except due to mypy false positive about redefining the function
    if sys.version_info >= (3, 12):
        from itertools import batched  # new in Python 3.12
    else:  # pragma: nocover
        from collections.abc import Iterator
        from itertools import islice

        def batched(iterable: Iterator[tuple], n: int) -> Iterator[tuple]:
            """Batch data from the iterable into tuples of length n.

            The last batch may be shorter than n.

            batched('ABCDEFG', 3) → ABC DEF G
            """
            if n < 1:
                msg = "n must be at least one"
                raise ValueError(msg)
            iterator = iter(iterable)
            while batch := tuple(islice(iterator, n)):
                yield batch

    uname = platform.uname()
    uname_system = uname.system
    uname_release = uname.release
    uname_machine = uname.machine

    configuration = run.configuration

    tool = tools.get_sourmash()
    _check_tool_version(logger, tool, configuration)

    config_id = configuration.configuration_id

    from pyani_plus.methods import sourmash  # lazy import

    sig_cache = cache / f"sourmash_k={configuration.kmersize}_{configuration.extra}"
    if not sig_cache.is_dir():
        msg = (
            f"Missing sourmash signatures directory '{sig_cache}'"
            f" - check cache setting '{cache}'."
        )
        log_sys_exit(logger, msg)

    try:
        db_entries: list[dict[str, str | float | None]] = []
        for batch in batched(
            sourmash.compute_sourmash_tile(
                logger,
                tool,
                {subject_hash} if subject_hash else set(query_hashes),
                set(query_hashes),
                sig_cache,
                tmp_dir,
            ),
            100000,
        ):
            logger.debug("Computed batch, about to log to database.")
            db_entries.extend(
                {
                    "query_hash": q,
                    "subject_hash": s,
                    "identity": max_containment,
                    "cov_query": q_containment,
                    "configuration_id": config_id,
                    "uname_system": uname_system,
                    "uname_release": uname_release,
                    "uname_machine": uname_machine,
                }
                for q, s, q_containment, max_containment in batch
            )
            export_json_db_entries(logger, json_filename, configuration, db_entries)
    except KeyboardInterrupt:  # pragma: no cover
        # Try to abort gracefully without wasting the work done.
        msg = f"Interrupted with {len(db_entries)} completed sourmash comparisons"
        logger.error(msg)  # noqa: TRY400
        run.status = "Worker interrupted"
        session.commit()

    try:
        export_json_db_entries(logger, json_filename, configuration, db_entries)
    except Exception:  # pragma: no cover
        logger.exception("Unexpected exception saving JSON:")
        return RECORDING_FAILED
    else:
        return 0


def compute_external_alignment(  # noqa: C901, PLR0912, PLR0913, PLR0915
    logger: logging.Logger,
    tmp_dir: Path,  # noqa: ARG001
    session: Session,
    run: db_orm.Run,
    json_filename: Path,
    fasta_dir: Path,  # noqa: ARG001
    hash_to_filename: dict[str, str],  # noqa: ARG001
    filename_to_hash: dict[str, str],  # noqa: ARG001
    query_hashes: dict[str, int],
    subject_hash: str,
    *,
    cache: Path = Path(),  # noqa: ARG001
) -> int:
    """Compute and log column of comparisons from given MSA to JSON.

    Will only look at query in query_hashes vs subject_hash, but will also
    record reciprocal comparison as this method is symmetric.
    """
    uname = platform.uname()
    uname_system = uname.system
    uname_release = uname.release
    uname_machine = uname.machine

    configuration = run.configuration
    config_id = configuration.configuration_id
    if configuration.method != "external-alignment":
        msg = f"Run-id {run.run_id} expected {configuration.method} results"
        log_sys_exit(logger, msg)
    if configuration.program:
        msg = f"configuration.program={configuration.program!r} unexpected"
        log_sys_exit(logger, msg)
    if configuration.version:
        msg = f"configuration.version={configuration.version!r} unexpected"
        log_sys_exit(logger, msg)
    if not configuration.extra:
        msg = "Missing configuration.extra setting"
        log_sys_exit(logger, msg)
    args = dict(_.split("=", 1) for _ in configuration.extra.split(";", 2))
    if list(args) != ["md5", "label", "alignment"]:
        msg = f"configuration.extra={configuration.extra!r} unexpected"
        log_sys_exit(logger, msg)

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
        msg = f"Treating {alignment} as relative to {url}"
        logger.debug(msg)
        alignment = Path(url[10:]).parent / alignment

    md5 = args["md5"]
    label = args["label"]
    del args

    from pyani_plus.methods.external_alignment import compute_external_alignment_column
    from pyani_plus.utils import file_md5sum, filename_stem

    msg = f"Parsing {alignment} (MD5={md5}, label={label})"
    logger.info(msg)
    if not alignment.is_file():
        msg = f"Missing alignment file {alignment}"
        log_sys_exit(logger, msg)
    if md5 != file_md5sum(alignment):
        msg = f"MD5 checksum of {alignment} didn't match."
        log_sys_exit(logger, msg)

    if label == "md5":
        mapping = lambda x: x  # noqa: E731
    elif label == "filename":
        mapping = {_.fasta_filename: _.genome_hash for _ in run.fasta_hashes}.get
    else:
        mapping = {
            filename_stem(_.fasta_filename): _.genome_hash for _ in run.fasta_hashes
        }.get

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
    try:
        db_entries = []
        for (
            q,
            s,
            identity,
            aln_length,
            sim_errors,
            cov_query,
            cov_subject,
        ) in compute_external_alignment_column(
            logger, subject_hash, set(query_hashes), alignment, mapping, label
        ):
            db_entries.append(
                {
                    "query_hash": q,
                    "subject_hash": s,
                    "identity": identity,
                    "aln_length": aln_length,
                    "sim_errors": sim_errors,
                    "cov_query": cov_query,
                    "cov_subject": cov_subject,
                    "configuration_id": config_id,
                    "uname_system": uname_system,
                    "uname_release": uname_release,
                    "uname_machine": uname_machine,
                }
            )
    except KeyboardInterrupt:  # pragma: no cover
        # Try to abort gracefully without wasting the work done.
        msg = f"Interrupted with {len(db_entries)} completed external-alignment comparisons"
        logger.error(msg)  # noqa: TRY400
        run.status = "Worker interrupted"
        session.commit()

    try:
        export_json_db_entries(logger, json_filename, configuration, db_entries)
    except Exception:  # pragma: no cover
        logger.exception("Unexpected exception saving JSON:")
        return RECORDING_FAILED
    else:
        return 0


if __name__ == "__main__":
    sys.exit(app())  # pragma: no cover
