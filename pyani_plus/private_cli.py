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

import platform
import sys
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
    method_dnadiff,
    method_fastani,
    method_sourmash,
)
from pyani_plus.public_cli import (
    OPT_ARG_TYPE_CREATE_DB,
    OPT_ARG_TYPE_FRAGSIZE,
    REQ_ARG_TYPE_DATABASE,
    REQ_ARG_TYPE_FASTA_DIR,
    REQ_ARG_TYPE_OUTDIR,
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
        help="Path to subject FASTA file",
        show_default=False,
        exists=True,
        dir_okay=False,
        file_okay=True,
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
    identity: Annotated[
        float,
        typer.Option(
            help="Percent identity",
            show_default=False,
            min=0.0,
            max=1.0,
        ),
    ],
    aln_length: Annotated[
        int, typer.Option(help="Alignment length", show_default=False)
    ],
    *,
    # Optional comparison table entries
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


# Ought we switch the command line arguments here to match fastANI naming?
# Note this omits mode
@app.command(rich_help_panel="Method specific logging")
def log_fastani(
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
    run_id: REQ_ARG_TYPE_RUN_ID,
    fastani: Annotated[
        Path,
        typer.Option(
            help="Path to fastANI output file",
            show_default=False,
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ],
    *,
    quiet: OPT_ARG_TYPE_QUIET = False,
) -> int:
    """Log fastANI pairwise comparison(s) to database.

    The associated configuration and genome entries must already exist.
    We expect this to be used on a row of the final matrix at a time,
    that is the output for many query genomes vs a single subject
    (reference) genome.
    """
    if database != ":memory:" and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    if not quiet:
        print(f"Logging fastANI comparison to {database}")
    session = db_orm.connect_to_db(database)
    run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()
    if run.configuration.method != "fastANI":
        msg = f"ERROR: Run-id {run_id} expected {run.configuration.method} results"
        sys.exit(msg)

    _check_tool_version(tools.get_fastani(), run.configuration)

    filename_to_hash = {_.fasta_filename: _.genome_hash for _ in run.fasta_hashes}

    for (
        query_hash,
        subject_hash,
        identity,
        orthologous_matches,
        fragments,
    ) in method_fastani.parse_fastani_file(fastani, filename_to_hash):
        # Allow for variation in the folder part of the filenames (e.g. relative paths)
        db_orm.db_comparison(
            session,
            configuration_id=run.configuration_id,
            query_hash=query_hash,
            subject_hash=subject_hash,
            identity=identity,
            aln_length=round(
                run.configuration.fragsize * orthologous_matches
            ),  # proxy value,
            sim_errors=fragments - orthologous_matches,  # proxy value, not bp,
            cov_query=float(orthologous_matches) / fragments,  # an approximation,
            cov_subject=None,
        )

    session.commit()
    return 0


@app.command()
def fragment_fasta(
    fasta: REQ_ARG_TYPE_FASTA_FILES,
    outdir: REQ_ARG_TYPE_OUTDIR,
    *,
    fragsize: OPT_ARG_TYPE_FRAGSIZE = method_anib.FRAGSIZE,
    quiet: OPT_ARG_TYPE_QUIET = False,
) -> int:
    """Fragment FASTA files into subsequences of up to the given size.

    The output files are named ``<stem>-fragmented.fna`` regardless of the
    input file extension (typically ``.fna``, ``.fa`` or ``.fasta``). If
    they already exist, they will be overwritten.
    """
    if not outdir.is_dir():
        sys.exit(f"ERROR: outdir {outdir} should be a directory")
    fragmented_files = method_anib.fragment_fasta_files(fasta, outdir, fragsize)
    if not quiet:
        print(f"Fragmented {len(fragmented_files)} files")
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
        cov_query=float(query_aligned_bases) / query.length,
        cov_subject=float(subject_aligned_bases) / subject.length,
    )

    session.commit()
    return 0


# Note this omits kmersize, minmatch, mode
@app.command(rich_help_panel="Method specific logging")
def log_anib(  # noqa: PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
    run_id: REQ_ARG_TYPE_RUN_ID,
    query_fasta: REQ_ARG_TYPE_QUERY_FASTA,
    subject_fasta: REQ_ARG_TYPE_SUBJECT_FASTA,
    blastn: Annotated[
        Path,
        typer.Option(
            help="Path to blastn TSV output file",
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ],
    *,
    quiet: OPT_ARG_TYPE_QUIET = False,
) -> int:
    """Log single ANIb pairwise comparison (with blastn) to database.

    The associated configuration and genome entries must already exist.
    """
    used_query, used_subject = blastn.stem.split("_vs_")
    if used_query != query_fasta.stem:
        sys.exit(
            f"ERROR: Given --query-fasta {query_fasta} but query in blastn filename was {used_query}"
        )
    if used_subject != subject_fasta.stem:
        sys.exit(
            f"ERROR: Given --subject-fasta {subject_fasta} but subject in blastn filename was {used_subject}"
        )

    if database != ":memory:" and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    if not quiet:
        print(f"Logging ANIb comparison to {database}")
    session = db_orm.connect_to_db(database)
    run, query_md5, subject_md5 = _lookup_run_query_subject(
        session, run_id, query_fasta, subject_fasta
    )
    if run.configuration.method != "ANIb":
        msg = f"ERROR: Run-id {run_id} expected {run.configuration.method} results"
        sys.exit(msg)

    _check_tool_version(tools.get_blastn(), run.configuration)

    # Need genome lengths for coverage (fragmented FASTA irrelevant):
    query = db_orm.db_genome(session, query_fasta, query_md5, create=False)
    subject = db_orm.db_genome(session, subject_fasta, subject_md5, create=False)

    identity, aln_length, sim_errors = method_anib.parse_blastn_file(blastn)

    db_orm.db_comparison(
        session,
        configuration_id=run.configuration_id,
        query_hash=query_md5,
        subject_hash=subject_md5,
        identity=identity,
        aln_length=aln_length,
        sim_errors=sim_errors,
        cov_query=float(aln_length) / query.length,
        cov_subject=float(aln_length) / subject.length,
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
        aln_length=aligned_bases_with_gaps - gap_lengths,
        sim_errors=round((aligned_bases_with_gaps - gap_lengths) * (1 - identity)),
        cov_query=(aligned_bases_with_gaps - gap_lengths) / query.length,
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


if __name__ == "__main__":
    sys.exit(app())  # pragma: no cover
