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
"""Implements the private command line interface (CLI) used internally.

The commands defined here are intended to be used from within pyANI-plus via
snakemake, for example from worker nodes, to log results to the database.
"""

import os
import platform
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Annotated

import typer
from rich.progress import Progress
from sqlalchemy.dialects.sqlite import insert as sqlite_insert
from sqlalchemy.orm import Session

from pyani_plus import PROGRESS_BAR_COLUMNS, db_orm, tools
from pyani_plus.methods import (
    method_anib,
    method_anim,
    method_branchwater,
    method_dnadiff,
    method_fastani,
    method_sourmash,
)
from pyani_plus.public_cli import (
    OPT_ARG_TYPE_CREATE_DB,
    OPT_ARG_TYPE_TEMP,
    REQ_ARG_TYPE_DATABASE,
    REQ_ARG_TYPE_FASTA_DIR,
)
from pyani_plus.utils import check_fasta, file_md5sum

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


def _lookup_run_query_subject(
    session: Session, run_id: int, query_fasta: Path, subject_fasta: Path
) -> tuple[db_orm.Run, str, str]:
    """Find run row in ORM, and lookup MD5 of query and subject FASTA files."""
    run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()
    query_md5 = (
        session.query(db_orm.RunGenomeAssociation)
        .where(db_orm.RunGenomeAssociation.run_id == run_id)
        .where(db_orm.RunGenomeAssociation.fasta_filename == query_fasta.name)
        .one()
        .genome_hash
    )
    subject_md5 = (
        session.query(db_orm.RunGenomeAssociation)
        .where(db_orm.RunGenomeAssociation.run_id == run_id)
        .where(db_orm.RunGenomeAssociation.fasta_filename == subject_fasta.name)
        .one()
        .genome_hash
    )
    return run, query_md5, subject_md5


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
    query_hashes: set[str] = set(hash_to_filename).difference(
        comp.query_hash
        for comp in run.comparisons().where(
            db_orm.Comparison.subject_hash == subject_hash
        )
    )

    if not query_hashes:
        if not quiet:
            print(f"INFO: No {method} comparisons needed against {subject_hash}")
        return 0

    # Will probably want to move each of these functions to the relevant method module...
    try:
        compute = {
            "fastANI": fastani,
            "ANIb": anib,
        }[method]
    except KeyError:
        msg = f"ERROR: Unknown method {method} for run-id {run_id} in {database}"
        sys.exit(msg)

    if temp:
        # Use the specified temp-directory (and do not clean up)
        return compute(
            temp,
            session,
            run,
            hash_to_filename,
            filename_to_hash,
            query_hashes,
            subject_hash,
            quiet=quiet,
        )
    # Use a system temp-directory (and do clean up)
    with tempfile.TemporaryDirectory() as sys_temp:
        return compute(
            Path(sys_temp),
            session,
            run,
            hash_to_filename,
            filename_to_hash,
            query_hashes,
            subject_hash,
            quiet=quiet,
        )


def fastani(  # noqa: PLR0913
    tmp_dir: Path,
    session: Session,
    run: db_orm.Run,
    hash_to_filename: dict[str, Path],
    filename_to_hash: dict[str, str],
    query_hashes: set[str],
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

    fasta_dir = Path(run.fasta_directory)
    tmp_output = tmp_dir / f"queries_vs_{subject_hash}.csv"
    tmp_queries = tmp_dir / f"queries_vs_{subject_hash}.txt"
    with tmp_queries.open("w") as handle:
        for query_hash in query_hashes:
            handle.write(f"{fasta_dir / hash_to_filename[query_hash]}\n")

    if not quiet:
        print(
            f"INFO: Calling fastANI for {len(query_hashes)} queries vs {subject_hash}"
        )

    # This will hide the stdout/stderr, unless it fails in which case we
    # get to see some of it via the default exception message:
    subprocess.check_output(
        [
            tool.exe_path,
            "--ql",
            str(tmp_queries),
            "-r",
            str(fasta_dir / hash_to_filename[subject_hash]),
            # Send to file or just capture stdout?
            "-o",
            str(tmp_output),
            "--fragLen",
            str(run.configuration.fragsize),
            "-k",
            str(run.configuration.kmersize),
            "--minFraction",
            str(run.configuration.minmatch),
        ],
        stderr=subprocess.STDOUT,
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
                "aln_length": round(run.configuration.fragsize * orthologous_matches),
                "sim_errors": fragments - orthologous_matches,
                "cov_query": float(orthologous_matches) / fragments,
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
            ) in method_fastani.parse_fastani_file(tmp_output, filename_to_hash)
        ],
    )
    session.commit()
    return 0


@app.command(rich_help_panel="Method specific logging")
def log_anim(  # noqa: PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
    run_id: REQ_ARG_TYPE_RUN_ID,
    query_fasta: REQ_ARG_TYPE_QUERY_FASTA,
    subject_fasta: REQ_ARG_TYPE_SUBJECT_FASTA,
    deltafilter: Annotated[
        Path,
        typer.Option(
            help="Path to deltafilter output file",
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ],
    *,
    quiet: OPT_ARG_TYPE_QUIET = False,
    # Don't use any of fragsize, maxmatch, kmersize, minmatch (configuration table entries)
) -> int:
    """Log single ANIm pairwise comparison (with nucmer) to database.

    The associated configuration and genome entries must already exist.
    """
    used_query, used_subject = deltafilter.stem.split("_vs_")
    if used_query != query_fasta.stem:
        sys.exit(
            f"ERROR: Given --query-fasta {query_fasta} but query in deltafilter filename was {used_query}"
        )
    if used_subject != subject_fasta.stem:
        sys.exit(
            f"ERROR: Given --subject-fasta {subject_fasta} but subject in deltafilter filename was {used_subject}"
        )

    if database != ":memory:" and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    if not quiet:
        print(f"Logging ANIm to {database}")
    session = db_orm.connect_to_db(database)
    run, query_md5, subject_md5 = _lookup_run_query_subject(
        session, run_id, query_fasta, subject_fasta
    )
    if run.configuration.method != "ANIm":
        msg = f"ERROR: Run-id {run_id} expected {run.configuration.method} results"
        sys.exit(msg)

    _check_tool_version(tools.get_nucmer(), run.configuration)

    # Need genome lengths for coverage:
    query = db_orm.db_genome(session, query_fasta, query_md5, create=False)
    subject = db_orm.db_genome(session, subject_fasta, subject_md5, create=False)

    query_aligned_bases, subject_aligned_bases, identity, sim_errors = (
        method_anim.parse_delta(deltafilter)
    )

    db_orm.db_comparison(
        session,
        configuration_id=run.configuration_id,
        query_hash=query_md5,
        subject_hash=subject_md5,
        identity=identity,
        aln_length=query_aligned_bases,
        sim_errors=sim_errors,
        cov_query=None
        if query_aligned_bases is None
        else float(query_aligned_bases) / query.length,
        cov_subject=None
        if subject_aligned_bases is None
        else float(subject_aligned_bases) / subject.length,
    )

    session.commit()
    return 0


def anib(  # noqa: PLR0913
    tmp_dir: Path,
    session: Session,
    run: db_orm.Run,
    hash_to_filename: dict[str, Path],
    filename_to_hash: dict[str, str],  # noqa: ARG001
    query_hashes: set[str],
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
    subject_length = (
        session.query(db_orm.Genome)
        .where(db_orm.Genome.genome_hash == subject_hash)
        .one()
        .length
    )
    subject_stem = Path(hash_to_filename[subject_hash]).stem
    outfmt = "6 " + " ".join(method_anib.BLAST_COLUMNS)

    fasta_dir = Path(run.fasta_directory)
    tmp_db = tmp_dir / f"{subject_stem}"  # prefix for BLAST DB

    if not quiet:
        print(f"INFO: Calling makeblastdb for {subject_stem}")

    # This will hide the stdout/stderr, unless it fails in which case we
    # get to see some of it via the default exception message:
    subprocess.check_output(
        [
            str(tools.get_makeblastdb().exe_path),
            "-in",
            str(fasta_dir / hash_to_filename[subject_hash]),
            "-input_type",
            "fasta",
            "-dbtype",
            "nucl",
            "-title",
            subject_hash,
            "-out",
            tmp_db,
        ],
        stderr=subprocess.STDOUT,
    )

    for query_hash in query_hashes:
        query_stem = Path(hash_to_filename[query_hash]).stem
        tmp_tsv = tmp_dir / f"{query_stem}_vs_{subject_stem}.tsv"

        # Potential race condition if other columns are being computed with the
        # same tmp_dir - so give the fragments file a unique name using PID:
        tmp_frag_query = (
            tmp_dir / f"{query_stem}-fragments-{fragsize}-pid{os.getpid()}.fna"
        )

        method_anib.fragment_fasta_file(
            fasta_dir / hash_to_filename[query_hash],
            tmp_frag_query,
            fragsize,
        )

        if not quiet:
            print(f"INFO: Calling blastn for {query_stem} vs {subject_stem}")

        # This will hide the stdout/stderr, unless it fails in which case we
        # get to see some of it via the default exception message:
        subprocess.check_output(
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
            stderr=subprocess.STDOUT,
        )

        identity, aln_length, sim_errors = method_anib.parse_blastn_file(tmp_tsv)

        query_length = (
            session.query(db_orm.Genome)
            .where(db_orm.Genome.genome_hash == query_hash)
            .one()
            .length
        )

        session.execute(
            sqlite_insert(db_orm.Comparison).on_conflict_do_nothing(),
            [
                {
                    "query_hash": query_hash,
                    "subject_hash": subject_hash,
                    "identity": identity,
                    "aln_length": aln_length,
                    "sim_errors": sim_errors,
                    "cov_query": float(aln_length) / query_length,
                    "cov_subject": float(aln_length) / subject_length,
                    "configuration_id": config_id,
                    "uname_system": uname_system,
                    "uname_release": uname_release,
                    "uname_machine": uname_machine,
                }
            ],
        )

        session.commit()
    return 0


@app.command(rich_help_panel="Method specific logging")
def log_dnadiff(  # noqa: PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
    run_id: REQ_ARG_TYPE_RUN_ID,
    query_fasta: REQ_ARG_TYPE_QUERY_FASTA,
    subject_fasta: REQ_ARG_TYPE_SUBJECT_FASTA,
    mcoords: Annotated[
        Path,
        typer.Option(
            help="Path to show-coords (.mcoords) output file",
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ],
    qdiff: Annotated[
        Path,
        typer.Option(
            help="Path to show-diff (.qdiff) output file",
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ],
    *,
    quiet: OPT_ARG_TYPE_QUIET = False,
) -> int:
    """Log single dnadiff pairwise comparison (with nucmer) to database.

    The associated configuration and genome entries must already exist.
    """
    # As with other methods, we need to verify that the provided query/subject sequences
    # match those used to generate the mcoords and qdiff files.
    # Checking if the query and subject used to generate mcoords files match those for
    # qdiff files might be unnecessary, as an incorrect query or subject will raise an error in lines 500-618.
    used_query_mcoords, used_subject_mcoords = mcoords.stem.split("_vs_")
    if used_query_mcoords != query_fasta.stem:
        sys.exit(
            f"ERROR: Given --query-fasta {query_fasta} but query in mcoords filename was {used_query_mcoords}"
        )
    if used_subject_mcoords != subject_fasta.stem:
        sys.exit(
            f"ERROR: Given --subject-fasta {subject_fasta} but subject in mcoords filename was {used_subject_mcoords}"
        )

    used_query_qdiff, used_subject_qdiff = qdiff.stem.split("_vs_")
    if used_query_qdiff != query_fasta.stem:
        sys.exit(
            f"ERROR: Given --query-fasta {query_fasta} but query in qdiff filename was {used_query_qdiff}"
        )
    if used_subject_qdiff != subject_fasta.stem:
        sys.exit(
            f"ERROR: Given --subject-fasta {subject_fasta} but subject in qdiff filename was {used_subject_qdiff}"
        )

    if database != ":memory:" and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    if not quiet:
        print(f"Logging dnadiff to {database}")
    session = db_orm.connect_to_db(database)
    run, query_md5, subject_md5 = _lookup_run_query_subject(
        session, run_id, query_fasta, subject_fasta
    )
    if run.configuration.method != "dnadiff":
        msg = f"ERROR: Run-id {run_id} expected {run.configuration.method} results"
        sys.exit(msg)

    _check_tool_version(tools.get_nucmer(), run.configuration)

    query = db_orm.db_genome(session, query_fasta, query_md5, create=False)

    identity, aligned_bases_with_gaps = method_dnadiff.parse_mcoords(mcoords)
    gap_lengths = method_dnadiff.parse_qdiff(qdiff)

    db_orm.db_comparison(
        session,
        configuration_id=run.configuration_id,
        query_hash=query_md5,
        subject_hash=subject_md5,
        identity=identity,
        # For comparisons of closely related genomes, qdiff files might
        # be empty as there are no gaps in the alignments. In this case, we
        # want to treat gap_lengths as 0. In cases of comparisons
        # of distantly related genomes, we report gap_lengths as None.
        aln_length=(
            None
            if gap_lengths is None and aligned_bases_with_gaps is None
            else (
                (aligned_bases_with_gaps or 0)
                - (gap_lengths if gap_lengths is not None else 0)
            )
        ),
        sim_errors=(
            None
            if identity is None or aligned_bases_with_gaps is None
            else round(
                (
                    (aligned_bases_with_gaps or 0)
                    - (gap_lengths if gap_lengths is not None else 0)
                )
                * (1 - identity)
            )
        ),
        cov_query=(
            None
            if aligned_bases_with_gaps is None or query.length == 0
            else (
                (aligned_bases_with_gaps or 0)
                - (gap_lengths if gap_lengths is not None else 0)
            )
            / query.length
        ),
        cov_subject=None,  # Leaving this as None for now (need rdiff files to calculate this)
    )

    session.commit()
    return 0


@app.command(rich_help_panel="Method specific logging")
def log_sourmash(
    database: REQ_ARG_TYPE_DATABASE,
    run_id: REQ_ARG_TYPE_RUN_ID,
    compare: Annotated[
        Path,
        typer.Option(
            help="Sourmash compare all-vs-all CSV output file",
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
                "identity": identity,
                "configuration_id": config_id,
                "uname_system": uname_system,
                "uname_release": uname_release,
                "uname_machine": uname_machine,
            }
            for (
                query_hash,
                subject_hash,
                identity,
            ) in method_sourmash.parse_sourmash_compare_csv(compare, filename_to_hash)
        ],
    )

    session.commit()
    return 0


@app.command(rich_help_panel="Method specific logging")
def log_branchwater(
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
        print(f"Logging branchwater to {database}")
    session = db_orm.connect_to_db(database)
    run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()
    if run.configuration.method != "branchwater":
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
                "identity": identity,
                "configuration_id": config_id,
                "uname_system": uname_system,
                "uname_release": uname_release,
                "uname_machine": uname_machine,
            }
            for (
                query_hash,
                subject_hash,
                identity,
            ) in method_branchwater.parse_sourmash_manysearch_csv(
                manysearch, filename_to_hash
            )
        ],
    )

    session.commit()
    return 0


if __name__ == "__main__":
    sys.exit(app())  # pragma: no cover
