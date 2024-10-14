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

import sys
from pathlib import Path
from typing import Annotated

import typer
from rich.progress import track

from pyani_plus import db_orm, tools
from pyani_plus.methods import method_anib, method_fastani
from pyani_plus.utils import file_md5sum

app = typer.Typer()


@app.command()
def log_configuration(  # noqa: PLR0913
    database: Annotated[
        str,
        typer.Option(help="Path to pyANI-plus SQLite3 database", show_default=False),
    ],
    method: Annotated[str, typer.Option(help="Method, e.g. ANIm", show_default=False)],
    program: Annotated[
        str, typer.Option(help="Program, e.g. nucmer", show_default=False)
    ],
    version: Annotated[
        str, typer.Option(help="Program version, e.g. 3.1", show_default=False)
    ],
    fragsize: Annotated[
        int | None, typer.Option(help="Optional method fragment size")
    ] = None,
    maxmatch: Annotated[
        bool | None, typer.Option(help="Comparison method max-match")
    ] = None,
    kmersize: Annotated[
        int | None, typer.Option(help="Comparison method k-mer size")
    ] = None,
    minmatch: Annotated[
        float | None, typer.Option(help="Comparison method min-match")
    ] = None,
    create_db: Annotated[  # noqa: FBT002
        bool, typer.Option(help="Create database if does not exist")
    ] = False,
) -> int:
    """Log this as a known configuration in the database.

    Any pre-existing configuration entry is left as is.
    """
    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist, but not using --create-db"
        sys.exit(msg)

    print(f"Logging to {database}")  # noqa: T201
    session = db_orm.connect_to_db(database)
    config = db_orm.add_configuration(
        session, method, program, version, fragsize, maxmatch, kmersize, minmatch
    )
    session.commit()  # should be redundant
    print(  # noqa: T201
        f"Configuration identifier {config.configuration_id}"
    )
    session.close()

    return 0


@app.command()
def log_genome(
    fasta: Annotated[
        list[str], typer.Argument(help="Path(s) to FASTA file(s)", show_default=False)
    ],
    database: Annotated[
        str,
        typer.Option(help="Path to pyANI-plus SQLite3 database", show_default=False),
    ],
    create_db: Annotated[  # noqa: FBT002
        bool, typer.Option(help="Create database if does not exist")
    ] = False,
) -> int:
    """For given FASTA file(s), compute their MD5 checksum, and log them in the database.

    Any pre-existing duplicate FASTA entries are left as is.
    """
    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist, but not using --create-db"
        sys.exit(msg)

    print(f"Logging to {database}")  # noqa: T201
    session = db_orm.connect_to_db(database)

    file_total = 0
    if fasta:
        for filename in track(fasta, description="Processing..."):
            file_total += 1
            md5 = file_md5sum(filename)
            db_orm.add_genome(session, filename, md5)
    session.commit()
    session.close()
    print(  # noqa: T201
        f"Processed {file_total} FASTA files"
    )

    return 0


@app.command()
def log_run(  # noqa: PLR0913
    fasta: Annotated[
        list[str], typer.Argument(help="Path(s) to FASTA file(s)", show_default=False)
    ],
    database: Annotated[
        str,
        typer.Option(help="Path to pyANI-plus SQLite3 database", show_default=False),
    ],
    # These are for the run table:
    cmdline: Annotated[str, typer.Option(help="Run command line", show_default=False)],
    status: Annotated[str, typer.Option(help="Run status", show_default=False)],
    name: Annotated[str, typer.Option(help="Run name", show_default=False)],
    # These are all for the configuration table:
    method: Annotated[str, typer.Option(help="Comparison method", show_default=False)],
    program: Annotated[
        str, typer.Option(help="Comparison program name", show_default=False)
    ],
    version: Annotated[
        str, typer.Option(help="Comparison program version", show_default=False)
    ],
    fragsize: Annotated[
        int | None, typer.Option(help="Comparison method fragment size")
    ] = None,
    maxmatch: Annotated[
        bool | None, typer.Option(help="Comparison method max-match")
    ] = None,
    kmersize: Annotated[
        int | None, typer.Option(help="Comparison method k-mer size")
    ] = None,
    minmatch: Annotated[
        float | None, typer.Option(help="Comparison method min-match")
    ] = None,
    create_db: Annotated[  # noqa: FBT002
        bool, typer.Option(help="Create database if does not exist")
    ] = False,
) -> int:
    """Log a run (and if need be, associated configuration, genome and comparison rows).

    There is currently no easy way to update an existing run (e.g. once more
    comparisons have been completed and you want to refresh the cached matrices
    and update the run status).
    """
    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist, but not using --create-db"
        sys.exit(msg)

    print(f"Logging to {database}")  # noqa: T201
    session = db_orm.connect_to_db(database)

    # Reuse existing config, or log a new one
    config = db_orm.add_configuration(
        session, method, program, version, fragsize, maxmatch, kmersize, minmatch
    )

    genomes = []
    if fasta:
        # Reuse existing genome entries and/or log new ones
        for filename in track(fasta, description="Processing..."):
            md5 = file_md5sum(filename)
            genomes.append(db_orm.add_genome(session, filename, md5))

    run = db_orm.add_run(
        session, config, cmdline, status, name, date=None, genomes=genomes
    )
    run.cache_comparisons()
    run_id = run.run_id

    session.commit()
    session.close()
    print(  # noqa: T201
        f"Run identifier {run_id}"
    )

    return 0


@app.command()
def log_comparison(  # noqa: PLR0913
    database: Annotated[
        str,
        typer.Option(help="Path to pyANI-plus SQLite3 database", show_default=False),
    ],
    # These are for the comparison table
    query_fasta: Annotated[
        Path, typer.Option(help="Path to query FASTA file", show_default=False)
    ],
    subject_fasta: Annotated[
        Path, typer.Option(help="Path to subject FASTA file", show_default=False)
    ],
    identity: Annotated[
        float,
        typer.Option(help="Percent identity (float from 0 to 1)", show_default=False),
    ],
    aln_length: Annotated[
        int, typer.Option(help="Alignment length", show_default=False)
    ],
    # These are all for the configuration table:
    method: Annotated[str, typer.Option(help="Comparison method", show_default=False)],
    program: Annotated[
        str, typer.Option(help="Comparison program name", show_default=False)
    ],
    version: Annotated[
        str, typer.Option(help="Comparison program version", show_default=False)
    ],
    fragsize: Annotated[
        int | None, typer.Option(help="Comparison method fragment size")
    ] = None,
    maxmatch: Annotated[
        bool | None, typer.Option(help="Comparison method max-match")
    ] = None,
    kmersize: Annotated[
        int | None, typer.Option(help="Comparison method k-mer size")
    ] = None,
    minmatch: Annotated[
        float | None, typer.Option(help="Comparison method min-match")
    ] = None,
    # Optional comparison table entries
    sim_errors: Annotated[int | None, typer.Option(help="Alignment length")] = None,
    cov_query: Annotated[float | None, typer.Option(help="Alignment length")] = None,
    cov_subject: Annotated[float | None, typer.Option(help="Alignment length")] = None,
    create_db: Annotated[  # noqa: FBT002
        bool, typer.Option(help="Create database if does not exist")
    ] = False,
) -> int:
    """Log a single pyANI-plus pairwise comparison to the database."""
    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist, but not using --create-db"
        sys.exit(msg)

    print(f"Logging to {database}")  # noqa: T201
    session = db_orm.connect_to_db(database)

    config = db_orm.add_configuration(
        session, method, program, version, fragsize, maxmatch, kmersize, minmatch
    )

    query_md5 = file_md5sum(query_fasta)
    db_orm.add_genome(session, query_fasta, query_md5)

    subject_md5 = file_md5sum(subject_fasta)
    db_orm.add_genome(session, subject_fasta, subject_md5)

    db_orm.add_comparison(
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
# Note this omits maxmatch
@app.command()
def log_fastani(  # noqa: PLR0913
    database: Annotated[
        str,
        typer.Option(help="Path to pyANI-plus SQLite3 database", show_default=False),
    ],
    # These are for the comparison table
    query_fasta: Annotated[
        Path, typer.Option(help="Path to query FASTA file", show_default=False)
    ],
    subject_fasta: Annotated[
        Path, typer.Option(help="Path to subject FASTA file", show_default=False)
    ],
    fastani: Annotated[
        Path, typer.Option(help="Path to fastANI output file", show_default=False)
    ],
    # These are all for the configuration table:
    fragsize: Annotated[
        int, typer.Option(help="Comparison method fragment size")
    ] = 3000,
    kmersize: Annotated[
        int | None, typer.Option(help="Comparison method k-mer size")
    ] = None,
    minmatch: Annotated[
        float | None, typer.Option(help="Comparison method min-match")
    ] = None,
    create_db: Annotated[  # noqa: FBT002
        bool, typer.Option(help="Create database if does not exist")
    ] = False,
) -> int:
    """Log a single pyANI-plus fastANI pairwise comparison to the database."""
    # Assuming this will match as expect this script to be called right
    # after the computation has finished (on the same machine)
    fastani_tool = tools.get_fastani()

    used_query, used_subject, identity, orthologous_matches, fragments = (
        method_fastani.parse_fastani_file(fastani)
    )
    # Allowing for some variation in the filename paths here... should we?
    if used_query.stem != query_fasta.stem:
        sys.exit(
            f"ERROR: Given --query-fasta {query_fasta} but query in fastANI file was {used_query}"
        )
    if used_subject.stem != subject_fasta.stem:
        sys.exit(
            f"ERROR: Given --subject-fasta {subject_fasta} but query in fastANI file was {used_subject}"
        )

    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist, but not using --create-db"
        sys.exit(msg)

    print(f"Logging to {database}")  # noqa: T201
    session = db_orm.connect_to_db(database)

    config = db_orm.add_configuration(
        session,
        method="fastANI",
        program=fastani_tool.exe_path.stem,
        version=fastani_tool.version,
        fragsize=fragsize,  # aka --fragLen
        maxmatch=None,
        kmersize=kmersize,  # aka --k
        minmatch=minmatch,  # aka --minFraction
    )

    query_md5 = file_md5sum(query_fasta)
    subject_md5 = file_md5sum(subject_fasta)

    # Should not be needed in standard workflow, but also ensures FASTA are in DB:
    db_orm.add_genome(session, query_fasta, query_md5)
    db_orm.add_genome(session, subject_fasta, subject_md5)

    estimated_cov_query = float(orthologous_matches) / fragments  # an approximation
    sim_errors = fragments - orthologous_matches  # proxy value, not bp
    estimated_aln_length = fragsize * orthologous_matches  # proxy value

    db_orm.add_comparison(
        session,
        configuration_id=config.configuration_id,
        query_hash=query_md5,
        subject_hash=subject_md5,
        identity=identity,
        aln_length=estimated_aln_length,
        sim_errors=sim_errors,
        cov_query=estimated_cov_query,
        cov_subject=None,
    )

    session.commit()
    return 0


@app.command()
def fragment_fasta(
    fasta: Annotated[
        list[Path],
        typer.Argument(help="Path to input FASTA file(s)", show_default=False),
    ],
    outdir: Annotated[Path, typer.Option(help="Output directory", show_default=False)],
    fragsize: Annotated[
        int, typer.Option(help="Fragment size (bp)")
    ] = method_anib.FRAGSIZE,
) -> int:
    """Fragment FASTA files into subsequences of up to the given size.

    The output files are named ``<stem>-fragmented.fna`` regardless of the
    input file extension (typically ``.fna``, ``.fa`` or ``.fasta``). If
    they already exist, they will be overwritten.
    """
    if not outdir.is_dir():
        sys.exit(f"ERROR: outdir {outdir} should be a directory")
    fragmented_files = method_anib.fragment_fasta_files(fasta, outdir, fragsize)
    print(f"Fragmented {len(fragmented_files)} files")  # noqa: T201
    return 0


# Note this omits kmersize, minmatch, maxmatch
@app.command()
def log_anib(  # noqa: PLR0913
    database: Annotated[str, typer.Option(help="Path to pyANI-plus SQLite3 database")],
    # These are for the comparison table
    query_fasta: Annotated[Path, typer.Option(help="Path to query FASTA file")],
    subject_fasta: Annotated[Path, typer.Option(help="Path to subject FASTA file")],
    blastn: Annotated[Path, typer.Option(help="Path to blastn TSV output file")],
    # These are all for the configuration table:
    fragsize: Annotated[
        int | None, typer.Option(help="Comparison method fragment size")
    ] = None,
    create_db: Annotated[  # noqa: FBT002
        bool, typer.Option(help="Create database if does not exist")
    ] = False,
) -> int:
    """Log a single pyANI-plus ANIb pairwise comparison (with blastn) to the database."""
    # Assuming this will match as expect this script to be called right
    # after the computation has finished (on the same machine)
    blastn_tool = tools.get_blastn()

    identity, aln_length, sim_errors = method_anib.parse_blastn_file(blastn)
    # Allowing for some variation in the filename paths here... should we?
    used_query, used_subject = blastn.stem.split("_vs_")
    if used_query != query_fasta.stem:
        sys.exit(
            f"ERROR: Given --query-fasta {query_fasta} but query in blastn filename was {used_query}"
        )
    if used_subject != subject_fasta.stem:
        sys.exit(
            f"ERROR: Given --subject-fasta {subject_fasta} but query in blastn filename was {used_subject}"
        )

    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist, but not using --create-db"
        sys.exit(msg)

    print(f"Logging to {database}")  # noqa: T201
    session = db_orm.connect_to_db(database)

    config = db_orm.add_configuration(
        session,
        method="ANIb",
        program=blastn_tool.exe_path.stem,
        version=blastn_tool.version,
        fragsize=fragsize,
        maxmatch=None,
        kmersize=None,
        minmatch=None,  # does a blastn filtering thresholds fit here?
    )

    # Need to lookup query length to compute query_cover (fragmented FASTA irrelevant:
    query_md5 = file_md5sum(query_fasta)
    query = db_orm.add_genome(session, query_fasta, query_md5)
    cov_query = float(aln_length) / query.length

    # Need to lookup subject length to compute subject_cover:
    subject_md5 = file_md5sum(subject_fasta)
    subject = db_orm.add_genome(session, subject_fasta, subject_md5)
    cov_subject = float(aln_length) / subject.length

    db_orm.add_comparison(
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


if __name__ == "__main__":
    sys.exit(app())  # pragma: no cover
