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

import subprocess
import sys
from pathlib import Path
from typing import Annotated

import typer
from rich.progress import Progress

from pyani_plus import PROGRESS_BAR_COLUMNS, db_orm, tools
from pyani_plus.methods import (
    method_anib,
    method_anim,
    method_dnadiff,
    method_fastani,
    method_sourmash,
)
from pyani_plus.public_cli import (
    OPT_ARG_TYPE_ANIM_MODE,
    OPT_ARG_TYPE_CREATE_DB,
    OPT_ARG_TYPE_FRAGSIZE,
    OPT_ARG_TYPE_KMERSIZE,
    OPT_ARG_TYPE_MINMATCH,
    OPT_ARG_TYPE_SOURMASH_MODE,
    REQ_ARG_TYPE_DATABASE,
    REQ_ARG_TYPE_FASTA_DIR,
    REQ_ARG_TYPE_OUTDIR,
    REQ_ARG_TYPE_RUN_NAME,
)
from pyani_plus.utils import check_fasta, file_md5sum

app = typer.Typer(
    context_settings={"help_option_names": ["-h", "--help"]},
)

OPT_ARG_TYPE_QUIET = Annotated[
    # Listing name(s) explicitly to avoid automatic matching --no-quiet
    bool, typer.Option("--quiet", help="Suppress any output except if fails")
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
    # These are all for the configuration table:
    method: REQ_ARG_TYPE_METHOD,
    program: REQ_ARG_TYPE_PROGRAM,
    version: REQ_ARG_TYPE_VERSION,
    *,
    fragsize: NONE_ARG_TYPE_FRAGSIZE = None,
    mode: NONE_ARG_TYPE_MODE = None,
    kmersize: NONE_ARG_TYPE_KMERSIZE = None,
    minmatch: NONE_ARG_TYPE_MINMATCH = None,
    extra: NONE_ARG_TYPE_EXTRA = None,
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
        create=False,
    )

    query_md5 = file_md5sum(query_fasta)
    db_orm.db_genome(session, query_fasta, query_md5)

    subject_md5 = file_md5sum(subject_fasta)
    db_orm.db_genome(session, subject_fasta, subject_md5)

    db_orm.db_comparison(
        session,
        configuration_id=config.configuration_id,
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


def _validate_method_inputs(  # noqa: PLR0913
    database: Path,
    method: str,
    query_fasta: Path,
    subject_fasta: Path,
    target_filenames: list[Path],
    tool: tools.ExternalToolData,
    *,
    fragsize: int | None = None,
    mode: str | None = None,
    kmersize: int | None = None,
    minmatch: float | None = None,
) -> tuple[db_orm.Session, db_orm.Configuration, str, str, bool]:
    """Validate method input arguments, compute query and subject hashes, and check if already in DB."""
    for filename in target_filenames:
        used_query, used_subject = filename.stem.split("_vs_")
        if used_query != query_fasta.stem:
            sys.exit(
                f"ERROR: Given --query-fasta {query_fasta} but query in target filename was {used_query}"
            )
        if used_subject != subject_fasta.stem:
            sys.exit(
                f"ERROR: Given --subject-fasta {subject_fasta} but subject in target filename was {used_subject}"
            )

    if database != ":memory:" and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    session = db_orm.connect_to_db(database)
    config = db_orm.db_configuration(
        session,
        method=method,
        program=tool.exe_path.stem,
        version=tool.version,
        fragsize=fragsize,
        mode=mode,
        kmersize=kmersize,
        minmatch=minmatch,
        create=False,
    )
    query_md5 = file_md5sum(query_fasta)
    subject_md5 = file_md5sum(subject_fasta)

    in_db = bool(
        session.query(db_orm.Comparison)
        .where(db_orm.Comparison.configuration_id == config.configuration_id)
        .where(db_orm.Comparison.query_hash == query_md5)
        .where(db_orm.Comparison.subject_hash == subject_md5)
        .count()
    )

    return session, config, query_md5, subject_md5, in_db


def _run_tool(args: list[str], *, quiet: bool) -> None:
    """Run a command line tool, hiding stdout/stderr unless if fails."""
    if not quiet:
        print("INFO: Running command: " + " ".join(args))
    try:
        subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as err:
        sys.stderr.write(f"ERROR: Command failed (return code {err.returncode}):\n")
        sys.stderr.write(" ".join(args) + "\n")
        sys.stderr.buffer.write(err.output)  # bytes!
        sys.exit(err.returncode)


def _run_tool_to_stdout(args: list[str], stdout_path: Path, *, quiet: bool) -> None:
    """Run a command line tool, stdout to given file, hiding stderr unless if fails."""
    if not quiet:
        print("INFO: Running command: " + " ".join(args) + f" > {stdout_path}")
    try:
        with stdout_path.open("w") as handle:
            subprocess.check_call(args, stdout=handle, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as err:
        sys.stderr.write(f"ERROR: Command failed (return code {err.returncode}):\n")
        sys.stderr.write(" ".join(args) + f" > {stdout_path}\n")
        sys.stderr.buffer.write(err.output)  # bytes!
        sys.exit(err.returncode)


# Ought we switch the command line arguments here to match fastANI naming?
# Note this omits mode
@app.command(rich_help_panel="ANI methods")
def fastani(  # noqa: PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
    query_fasta: REQ_ARG_TYPE_QUERY_FASTA,
    subject_fasta: REQ_ARG_TYPE_SUBJECT_FASTA,
    fastani: Annotated[
        Path,
        typer.Option(
            help="Path to fastANI output file (will be generated if missing)",
            show_default=False,
            dir_okay=False,
            file_okay=True,
        ),
    ],
    *,
    quiet: OPT_ARG_TYPE_QUIET = False,
    # These are all for the configuration table:
    fragsize: OPT_ARG_TYPE_FRAGSIZE = method_fastani.FRAG_LEN,
    kmersize: OPT_ARG_TYPE_KMERSIZE = method_fastani.KMER_SIZE,
    minmatch: OPT_ARG_TYPE_MINMATCH = method_fastani.MIN_FRACTION,
) -> int:
    """Compute a single fastANI pairwise comparison and log to database.

    The associated configuration and genome entries must already exist.
    If the fastANI output file already exists, it will be reused. If not,
    this will run the fastANI tool to generate that file.
    """
    tool = tools.get_fastani()
    session, config, query_md5, subject_md5, in_db = _validate_method_inputs(
        database,
        "fastANI",
        query_fasta,
        subject_fasta,
        [fastani],
        tool,
        fragsize=fragsize,
        kmersize=kmersize,
        minmatch=minmatch,
    )

    if in_db:
        if not quiet:
            print("Already in the database")
        session.close()
        return 0

    if not fastani.is_file():
        _run_tool(
            [
                str(tool.exe_path),
                "-q",
                str(query_fasta),
                "-r",
                str(subject_fasta),
                "-o",
                str(fastani),
                "--fragLen",
                str(fragsize),
                "-k",
                str(kmersize),
                "--minFraction",
                str(minmatch),
            ],
            quiet=quiet,
        )

    if not quiet:
        print(f"Parsing {fastani}")
    used_query_path, used_subject_path, identity, orthologous_matches, fragments = (
        method_fastani.parse_fastani_file(fastani)
    )
    # Allow for variation in the folder part of the filenames (e.g. relative paths)
    if used_query_path.stem != query_fasta.stem:
        sys.exit(
            f"ERROR: Given --query-fasta {query_fasta}"
            f" but query in fastANI file contents was {used_query_path}"
        )
    if used_subject_path.stem != subject_fasta.stem:
        sys.exit(
            f"ERROR: Given --subject-fasta {subject_fasta}"
            f" but subject in fastANI file contents was {used_subject_path}"
        )

    # We assume both genomes have been recorded, if not this will fail:
    db_orm.db_comparison(
        session,
        configuration_id=config.configuration_id,
        query_hash=query_md5,
        subject_hash=subject_md5,
        identity=identity,
        aln_length=round(fragsize * orthologous_matches),  # proxy value,
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


@app.command(rich_help_panel="ANI methods")
def anim(  # noqa: PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
    query_fasta: REQ_ARG_TYPE_QUERY_FASTA,
    subject_fasta: REQ_ARG_TYPE_SUBJECT_FASTA,
    delta: Annotated[
        Path,
        typer.Option(
            help="Path to nucmer delta output file",
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ],
    deltafilter: Annotated[
        Path,
        typer.Option(
            help="Path to deltafilter output file (will be generated if missing)",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    *,
    # Don't use any of fragsize, kmersize, minmatch (configuration table entries)
    mode: OPT_ARG_TYPE_ANIM_MODE = method_anim.MODE,
    quiet: OPT_ARG_TYPE_QUIET = False,
) -> int:
    """Compute a single ANIm pairwise comparison and log to database.

    The associated configuration and genome entries must already exist.
    """
    tool = tools.get_nucmer()
    tool_deltafilter = tools.get_delta_filter()
    if (
        delta.stem != deltafilter.stem
        or delta.suffix != ".delta"
        or deltafilter.suffix != ".filter"
    ):
        msg = "ERROR: Require the delta and delta-filter files to be named <stem>.delta and <stem>.filter"
        sys.exit(msg)
    session, config, query_md5, subject_md5, in_db = _validate_method_inputs(
        database,
        "ANIm",
        query_fasta,
        subject_fasta,
        [delta, deltafilter],
        tool,
        mode=mode.value,  # enum to string
    )

    if in_db:
        if not quiet:
            print("Already in the database")
        session.close()
        return 0

    if not deltafilter.is_file():
        _run_tool_to_stdout(
            [str(tool_deltafilter.exe_path), "-1", str(delta)],
            quiet=quiet,
            stdout_path=deltafilter,
        )

    if not quiet:
        print("Parsing ANIm comparison from nucmer via delta-filter")
    query_aligned_bases, subject_aligned_bases, identity, sim_errors = (
        method_anim.parse_delta(deltafilter)
    )

    # Need genome lengths for coverage
    query = db_orm.db_genome(session, query_fasta, query_md5, create=False)
    subject = db_orm.db_genome(session, subject_fasta, subject_md5, create=False)

    db_orm.db_comparison(
        session,
        configuration_id=config.configuration_id,
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
@app.command(rich_help_panel="ANI methods")
def anib(  # noqa: PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
    query_fasta: REQ_ARG_TYPE_QUERY_FASTA,
    subject_fasta: REQ_ARG_TYPE_SUBJECT_FASTA,
    fragments: Annotated[
        Path,
        typer.Option(
            help="Path to fragmented query file",
            show_default=False,
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ],
    blastdb: Annotated[
        Path,
        typer.Option(
            help="Path to subject BLAST database .njs file",
            show_default=False,
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ],
    blastn: Annotated[
        Path,
        typer.Option(
            help="Path to blastn TSV output file",
            show_default=False,
            dir_okay=False,
            file_okay=True,
        ),
    ],
    *,
    quiet: OPT_ARG_TYPE_QUIET = False,
    # These are all for the configuration table:
    fragsize: OPT_ARG_TYPE_FRAGSIZE = method_anib.FRAGSIZE,
) -> int:
    """Compute a single ANIb pairwise comparison and log to database.

    The associated configuration and genome entries must already exist.
    """
    tool = tools.get_blastn()
    if fragments.name != query_fasta.stem + "-fragments.fna":
        msg = (
            "ERROR - Fragmented version of query should be named"
            f" {query_fasta.stem}-fragments.fna, not {fragments.name}"
        )
        sys.exit(msg)
    if blastdb.name != subject_fasta.stem + ".njs":
        msg = f"ERROR - Require the subject BLASTDB to be named {subject_fasta.stem}.njs, not {blastdb.name}"
        sys.exit(msg)
    session, config, query_md5, subject_md5, in_db = _validate_method_inputs(
        database,
        "ANIb",
        query_fasta,
        subject_fasta,
        [blastn],
        tool,
        fragsize=fragsize,
    )

    if in_db:
        if not quiet:
            print("Already in the database")
        session.close()
        return 0

    if not blastn.is_file():
        _run_tool(
            [
                str(tool.exe_path),
                "-query",
                str(fragments),
                "-db",
                str(blastdb.with_suffix("")),  # drop the .njs
                "-out",
                str(blastn),
                "-task",
                "blastn",
                "-outfmt",
                "6 qseqid sseqid length mismatch pident nident qlen slen qstart qend sstart send positive ppos gaps",
                "-xdrop_gap_final",
                "150",
                "-dust",
                "no",
                "-evalue",
                "1e-15",
                "-max_target_seqs",
                "1",
            ],
            quiet=quiet,
        )

    if not quiet:
        print(f"Parsing ANIb comparison from {blastn}")
    identity, aln_length, sim_errors = method_anib.parse_blastn_file(blastn)

    # Need to lookup query and subject lengths for coverage (fragment length no relevant)
    query = db_orm.db_genome(session, query_fasta, query_md5, create=False)
    subject = db_orm.db_genome(session, subject_fasta, subject_md5, create=False)

    if not quiet:
        print(f"Logging ANIb comparison to {database}")
    db_orm.db_comparison(
        session,
        configuration_id=config.configuration_id,
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


@app.command(rich_help_panel="ANI methods")
def dnadiff(  # noqa: PLR0913, C901
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
    query_fasta: REQ_ARG_TYPE_QUERY_FASTA,
    subject_fasta: REQ_ARG_TYPE_SUBJECT_FASTA,
    deltafilter: Annotated[
        Path,
        typer.Option(
            help="Path to deltafilter (.filter) output file",
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ],
    mcoords: Annotated[
        Path,
        typer.Option(
            help="Path to show-coords (.mcoords) output file (will be generated if missing)",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    qdiff: Annotated[
        Path,
        typer.Option(
            help="Path to show-diff (.qdiff) output file (will be generated if missing)",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    *,
    quiet: OPT_ARG_TYPE_QUIET = False,
    # Should we add --maxmatch (nucmer) and -m (deltafilter) parameters?
    # These are default parameters used in workflows
) -> int:
    """Compute a single dnadiff pairwise comparison and log to database.

    The associated configuration and genome entries must already exist.
    """
    tool = tools.get_nucmer()
    tool_show_diff = tools.get_show_diff()
    tool_show_coords = tools.get_show_coords()
    if (
        not (deltafilter.stem == mcoords.stem == qdiff.stem)
        or deltafilter.suffix != ".filter"
        or mcoords.suffix != ".mcoords"
        or qdiff.suffix != ".qdiff"
    ):
        msg = (
            "ERROR: Require the delta-filter, qdiff, and mcoords files to be"
            "named <stem>.filter, <stem>.qdiff, and <stem>.mcoords"
        )
        sys.exit(msg)
    session, config, query_md5, subject_md5, in_db = _validate_method_inputs(
        database,
        "dnadiff",
        query_fasta,
        subject_fasta,
        [deltafilter, qdiff, mcoords],
        tool,
    )

    if in_db:
        if not quiet:
            print("Already in the database")
        session.close()
        return 0

    if not qdiff.is_file():
        _run_tool_to_stdout(
            [str(tool_show_diff.exe_path), "-qH", str(deltafilter)],
            quiet=quiet,
            stdout_path=qdiff,
        )
    if not mcoords.is_file():
        _run_tool_to_stdout(
            [str(tool_show_coords.exe_path), "-rclTH", str(deltafilter)],
            quiet=quiet,
            stdout_path=mcoords,
        )

    if not quiet:
        print("Parsing ANIm comparison from nucmer via delta-filter")
    identity, aligned_bases_with_gaps = method_dnadiff.parse_mcoords(mcoords)
    gap_lengths = method_dnadiff.parse_qdiff(qdiff)

    # As with other methods, we need to verify that the provided query/subject sequences
    # match those used to generate the mcoords and qdiff files.
    # Checking if the query and subject used to generate mcoords files match those for
    # qdiff files might be unnecessary, as an incorrect query or subject will raise an error in lines 500-618.
    used_query_mcoords, used_subject_mcoords = mcoords.stem.split("_vs_")
    if used_query_mcoords != query_fasta.stem:
        sys.exit(
            f"ERROR: Given --query-fasta {query_fasta} but query in target filename was {used_query_mcoords}"
        )
    if used_subject_mcoords != subject_fasta.stem:
        sys.exit(
            f"ERROR: Given --subject-fasta {subject_fasta} but subject in target filename was {used_subject_mcoords}"
        )

    used_query_qdiff, used_subject_qdiff = qdiff.stem.split("_vs_")
    if used_query_qdiff != query_fasta.stem:
        sys.exit(
            f"ERROR: Given --query-fasta {query_fasta} but query in target filename was {used_query_qdiff}"
        )
    if used_subject_qdiff != subject_fasta.stem:
        sys.exit(
            f"ERROR: Given --subject-fasta {subject_fasta} but subject in target filename was {used_subject_qdiff}"
        )

    if database != ":memory:" and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    if not quiet:
        print(f"Logging dnadiff to {database}")
    session = db_orm.connect_to_db(database)

    config = db_orm.db_configuration(
        session,
        method="dnadiff",
        program=tool.exe_path.stem,
        version=tool.version,
        create=False,
    )

    query_md5 = file_md5sum(query_fasta)
    subject_md5 = file_md5sum(subject_fasta)

    query = db_orm.db_genome(session, query_fasta, query_md5, create=False)

    db_orm.db_comparison(
        session,
        configuration_id=config.configuration_id,
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


@app.command(rich_help_panel="ANI methods")
def sourmash(  # noqa: PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
    query_fasta: REQ_ARG_TYPE_QUERY_FASTA,
    subject_fasta: REQ_ARG_TYPE_SUBJECT_FASTA,
    compare: Annotated[
        Path,
        typer.Option(
            help="Path to sourmash compare CSV output file",
            show_default=False,
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ],
    *,
    quiet: OPT_ARG_TYPE_QUIET = False,
    # # Don't use any of fragsize, minmatch (configuration table entries)
    mode: OPT_ARG_TYPE_SOURMASH_MODE = method_sourmash.MODE,
    kmersize: OPT_ARG_TYPE_KMERSIZE = method_sourmash.KMER_SIZE,
    extra: NONE_ARG_TYPE_EXTRA = f"scaled={method_sourmash.SCALED}",
    # Don't use any of fragsize, maxmatch, minmatch (configuration table entries)
) -> int:
    """Compute a single dnadiff pairwise comparison and log to database.

    The associated configuration and genome entries must already exist.
    """
    # Assuming this will match as expect this script to be called right
    # after the computation has finished (on the same machine)
    sourmash_tool = tools.get_sourmash()

    identity = method_sourmash.parse_compare(compare)

    used_query, used_subject = compare.stem.split("_vs_")
    if used_query != query_fasta.stem:
        sys.exit(
            f"ERROR: Given --query-fasta {query_fasta} but query in sourmash compare filename was {used_query}"
        )
    if used_subject != subject_fasta.stem:
        sys.exit(
            f"ERROR: Given --subject-fasta {subject_fasta} but subject in sourmash compare filename was {used_subject}"
        )

    if database != ":memory:" and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    if not quiet:
        print(f"Logging sourmash to {database}")
    session = db_orm.connect_to_db(database)

    config = db_orm.db_configuration(
        session,
        method="sourmash",
        program=sourmash_tool.exe_path.stem,
        version=sourmash_tool.version,
        mode=mode.value,  # turn enum into string
        kmersize=kmersize,
        extra=extra,
        create=False,
    )

    query_md5 = file_md5sum(query_fasta)
    subject_md5 = file_md5sum(subject_fasta)

    db_orm.db_comparison(
        session,
        configuration_id=config.configuration_id,
        query_hash=query_md5,
        subject_hash=subject_md5,
        identity=identity,
        # Question: We need to discuss how to calculate or report the below values.
        aln_length=None,
        sim_errors=None,
        cov_query=None,
        cov_subject=None,
    )

    session.commit()
    return 0


if __name__ == "__main__":
    sys.exit(app())  # pragma: no cover
