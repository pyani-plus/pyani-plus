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

import importlib
import platform
import sys
from itertools import product
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

REQ_ARG_TYPE_CACHE = Annotated[
    Path | None,
    typer.Option(
        help=(
            "Cache location (must be visible to cluster workers)."
            " Default .cache/<method>"
        ),
        metavar="DIRECTORY",
        # Not requiring this exists, not all methods use a cache
        dir_okay=True,
        file_okay=False,
    ),
]


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


@app.command(rich_help_panel="Main")
def prepare(
    database: REQ_ARG_TYPE_DATABASE,
    run_id: Annotated[
        int | None,
        typer.Option(help="Which run to prepare", show_default=False),
    ],
    cache: REQ_ARG_TYPE_CACHE = None,
    *,
    quiet: OPT_ARG_TYPE_QUIET = False,
) -> int:
    """Prepare any intermediate files prior to computing the ANI values.

    This requires you already have a run defined in the database, and would
    be followed by running the private CLI ``compute`` command.
    """
    # Should this be splittable for running on the cluster? I assume most
    # cases this is IO bound rather than CPU bound so is this helpful?
    if database != ":memory:" and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)
    session = db_orm.connect_to_db(database)
    run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()
    n = run.genomes.count()
    done = run.comparisons().count()
    if done == n**2:
        if not quiet:
            print(
                f"Skipping preparation, run already has all {n**2}={n}² pairwise values"
            )
        return 0
    config = run.configuration
    method = config.method
    module = importlib.import_module(f"pyani_plus.methods.method_{method.lower()}")
    if not hasattr(module, "prepare_genomes"):
        sys.stderr.write(f"No per-genome preparation required for {method}\n")
        return 0

    if cache is None:
        cache = Path(f".cache/{method}/")
        cache.mkdir(parents=True, exist_ok=True)

    with Progress(*PROGRESS_BAR_COLUMNS) as progress:
        for _ in progress.track(
            module.prepare_genomes(run, cache),
            description="Processing...  ",  # spaces to match "Indexing FASTAs" etc
            total=n,
        ):
            pass
    return 0


@app.command(rich_help_panel="Main")
def compute(  # noqa: C901, PLR0912, PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    run_id: Annotated[
        int | None,
        typer.Option(help="Which run to resume", show_default=False),
    ],
    cache: REQ_ARG_TYPE_CACHE = None,
    parts: Annotated[
        int,
        typer.Option(
            help=(
                "How many workers, or how many parts,"
                " is the computation to be split between?"
            ),
            min=1,
        ),
    ] = 1,
    task: Annotated[
        int,
        typer.Option(
            help=(
                "Given the computation is broken into parts, which one to compute?"
                " Should be between 0 and number of parts."
            ),
            min=0,
        ),
    ] = 0,
    *,
    quiet: OPT_ARG_TYPE_QUIET = False,
) -> int:
    """Compute the ANI values for a run already defined in the database.

    For some methods you must first have run the private CLI ``prepare``
    command, e.g. for ANIb this builds BLAST databases, or for sourmash
    this builds signature sketches.

    If you set parts equal to the number of genomes, this will make each
    row of the database into a task (i.e. all vs a single subject FASTA).
    """
    if not (0 <= task <= parts):
        msg = f"ERROR: Expect task in range 0 to {parts}, got {task}"
        sys.exit(msg)
    if task == parts:
        # This is to accommodate some flexibility in array job numbering
        task = 0
    if database != ":memory:" and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)
    uname = platform.uname()
    session = db_orm.connect_to_db(database)
    run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()
    n = run.genomes.count()
    done = run.comparisons().count()
    if not quiet:
        print(f"Run already has {done}/{n**2}={n}² pairwise values")
    if done == n**2:
        return 0
    fasta_dir = Path(run.fasta_directory)
    config = run.configuration
    method = config.method
    module = importlib.import_module(f"pyani_plus.methods.method_{method.lower()}")

    if cache is None:
        cache = Path(f".cache/{method}/")
    if hasattr(module, "prepare_genomes") and not cache.is_dir():
        msg = f"ERROR: {method} needs prepared files but cache {cache} directory does not exist"
        sys.exit(msg)

    hashes = {
        _.genome_hash: _.fasta_filename
        for _ in run.fasta_hashes.order_by(db_orm.RunGenomeAssociation.genome_hash)
    }

    for index, (query_md5, subject_md5) in enumerate(product(hashes, hashes)):
        if index % parts != task:
            if not quiet:
                print(f"Skipping {query_md5} vs {subject_md5} as not requested")
            continue
        if bool(
            session.query(db_orm.Comparison)
            .where(db_orm.Comparison.configuration_id == config.configuration_id)
            .where(db_orm.Comparison.query_hash == query_md5)
            .where(db_orm.Comparison.subject_hash == subject_md5)
            .count()
        ):
            if not quiet:
                print(f"Skipping {query_md5} vs {subject_md5} as already in DB")
            continue
        # Actually do the computation!
        if not quiet:
            print(f"Computing {query_md5} vs {subject_md5} now...")
        query = (
            session.query(db_orm.Genome)
            .where(db_orm.Genome.genome_hash == query_md5)
            .one()
        )
        subject = (
            session.query(db_orm.Genome)
            .where(db_orm.Genome.genome_hash == subject_md5)
            .one()
        )
        session.add(
            module.compute_pairwise_ani(
                uname,
                config,
                query_md5,
                fasta_dir / hashes[query_md5],
                query.length,
                subject_md5,
                fasta_dir / hashes[subject_md5],
                subject.length,
                cache,
            )
        )
        session.commit()
    return 0


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


# Ought we switch the command line arguments here to match fastANI naming?
# Note this omits mode
@app.command(rich_help_panel="Method specific logging")
def log_fastani(  # noqa: PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
    query_fasta: REQ_ARG_TYPE_QUERY_FASTA,
    subject_fasta: REQ_ARG_TYPE_SUBJECT_FASTA,
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
    # These are all for the configuration table:
    fragsize: OPT_ARG_TYPE_FRAGSIZE = method_fastani.FRAG_LEN,
    kmersize: OPT_ARG_TYPE_KMERSIZE = method_fastani.KMER_SIZE,
    minmatch: OPT_ARG_TYPE_MINMATCH = method_fastani.MIN_FRACTION,
) -> int:
    """Log single fastANI pairwise comparison to database.

    The associated configuration and genome entries must already exist.
    """
    # Assuming this will match as expect this script to be called right
    # after the computation has finished (on the same machine)
    fastani_tool = tools.get_fastani()

    used_query, used_subject = fastani.stem.split("_vs_")
    if used_query != query_fasta.stem:
        sys.exit(
            f"ERROR: Given --query-fasta {query_fasta} but query in fastANI filename was {used_query}"
        )
    if used_subject != subject_fasta.stem:
        sys.exit(
            f"ERROR: Given --subject-fasta {subject_fasta} but subject in fastANI filename was {used_subject}"
        )

    used_query_path, used_subject_path, identity, orthologous_matches, fragments = (
        method_fastani.parse_fastani_file(fastani)
    )
    # Allow for variation in the folder part of the filenames (e.g. relative paths)
    if used_query_path.stem != query_fasta.stem:
        sys.exit(
            f"ERROR: Given --query-fasta {query_fasta} but query in fastANI file contents was {used_query_path}"
        )
    if used_subject_path.stem != subject_fasta.stem:
        sys.exit(
            f"ERROR: Given --subject-fasta {subject_fasta} but subject in fastANI file contents was {used_subject_path}"
        )

    if database != ":memory:" and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist"
        sys.exit(msg)

    if not quiet:
        print(f"Logging fastANI comparison to {database}")
    session = db_orm.connect_to_db(database)

    config = db_orm.db_configuration(
        session,
        method="fastANI",
        program=fastani_tool.exe_path.stem,
        version=fastani_tool.version,
        fragsize=fragsize,  # aka --fragLen
        kmersize=kmersize,  # aka --k
        minmatch=minmatch,  # aka --minFraction
        create=False,
    )

    query_md5 = file_md5sum(query_fasta)
    subject_md5 = file_md5sum(subject_fasta)

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


@app.command(rich_help_panel="Method specific logging")
def log_anim(  # noqa: PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
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
    # Don't use any of fragsize, kmersize, minmatch (configuration table entries)
    mode: OPT_ARG_TYPE_ANIM_MODE = method_anim.MODE,
    quiet: OPT_ARG_TYPE_QUIET = False,
    # Don't use any of fragsize, maxmatch, kmersize, minmatch (configuration table entries)
) -> int:
    """Log single ANIm pairwise comparison (with nucmer) to database.

    The associated configuration and genome entries must already exist.
    """
    # Assuming this will match as expect this script to be called right
    # after the computation has finished (on the same machine)
    nucmer_tool = tools.get_nucmer()

    query_aligned_bases, subject_aligned_bases, identity, sim_errors = (
        method_anim.parse_delta(deltafilter)
    )

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

    config = db_orm.db_configuration(
        session,
        method="ANIm",
        program=nucmer_tool.exe_path.stem,
        version=nucmer_tool.version,
        mode=mode,
        create=False,
    )

    query_md5 = file_md5sum(query_fasta)
    subject_md5 = file_md5sum(subject_fasta)

    # Need genome lengths for coverage:
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
@app.command(rich_help_panel="Method specific logging")
def log_anib(  # noqa: PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
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
    # These are all for the configuration table:
    fragsize: OPT_ARG_TYPE_FRAGSIZE = method_anib.FRAGSIZE,
) -> int:
    """Log single ANIb pairwise comparison (with blastn) to database.

    The associated configuration and genome entries must already exist.
    """
    # Assuming this will match as expect this script to be called right
    # after the computation has finished (on the same machine)
    blastn_tool = tools.get_blastn()

    identity, aln_length, sim_errors = method_anib.parse_blastn_file(blastn)
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

    config = db_orm.db_configuration(
        session,
        method="ANIb",
        program=blastn_tool.exe_path.stem,
        version=blastn_tool.version,
        fragsize=fragsize,
        create=False,
    )

    query_md5 = file_md5sum(query_fasta)
    subject_md5 = file_md5sum(subject_fasta)

    # Need genome lengths for coverage (fragmented FASTA irrelevant):
    query = db_orm.db_genome(session, query_fasta, query_md5, create=False)
    subject = db_orm.db_genome(session, subject_fasta, subject_md5, create=False)

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


@app.command(rich_help_panel="Method specific logging")
def log_dnadiff(  # noqa: PLR0913
    database: REQ_ARG_TYPE_DATABASE,
    # These are for the comparison table
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
    # Should we add --maxmatch (nucmer) and -m (deltafilter) parameters?
    # These are default parameters used in workflows
) -> int:
    """Log single dnadiff pairwise comparison (with nucmer) to database.

    The associated configuration and genome entries must already exist.
    """
    # Assuming this will match as expect this script to be called right
    # after the computation has finished (on the same machine)
    # We don't actually call the tool dnadiff (which has its own version),
    # rather we call nucmer, delta-filter, show-diff and show-coords from mumer
    tool = tools.get_nucmer()

    identity, aligned_bases_with_gaps = method_dnadiff.parse_mcoords(mcoords)
    gap_lengths = method_dnadiff.parse_qdiff(qdiff)

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


@app.command(rich_help_panel="Method specific logging")
def log_sourmash(  # noqa: PLR0913
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
    """Log single sourmash pairwise comparison to database.

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
