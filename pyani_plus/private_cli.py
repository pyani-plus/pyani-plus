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
from pyani_plus.utils import file_md5sum

app = typer.Typer()


@app.command()
def log_genome(
    fasta: Annotated[list[str], typer.Argument(help="Path to FASTA file(s)")],
    database: Annotated[str, typer.Option(help="Path to pyANI-plus SQLite3 database")],
) -> int:
    """For given FASTA file(s), compute their MD5 checksum, and log them in the database.

    Any pre-existing duplicate FASTA entries are left as is.
    """
    print(f"Logging to {database}")  # noqa: T201
    db = Path(database)
    if not db.is_file():
        sys.exit(f"ERROR - Database {db} does not exist")
    session = db_orm.connect_to_db(db)

    file_total = 0
    file_skip = 0
    if fasta:
        for filename in track(fasta, description="Processing..."):
            file_total += 1
            md5 = file_md5sum(filename)
            if not db_orm.add_genome(session, filename, md5):
                file_skip += 1
    session.commit()
    session.close()
    print(  # noqa: T201
        f"Processed {file_total} FASTA files, skipped {file_skip}, recorded {file_total-file_skip}"
    )

    return 0


def _log_comparison(  # noqa: PLR0913
    session: db_orm.Session,
    config: db_orm.Configuration,
    query_fasta: Path,
    subject_fasta: Path,
    identity: float,
    aln_length: int,
    sim_errors: int | None = None,
) -> bool:
    """Record the two genomes and comparison itself.

    Will first compute the MD5 checksum for each FASTA file, and add them to the
    genomes table if required. Then using those and the given configuration ID,
    will add the comparison if required.

    This will call session.commit() unconditionally.

    Returns True for a novel comparison, False if already present.
    """
    query_md5 = file_md5sum(query_fasta)
    subject_md5 = file_md5sum(subject_fasta)

    if db_orm.add_genome(session, query_fasta, query_md5):
        print(f"Adding novel FASTA hash {query_md5} to database")  # noqa: T201
    if db_orm.add_genome(session, subject_fasta, subject_md5):
        print(f"Adding novel FASTA hash {subject_md5} to database")  # noqa: T201
    new = db_orm.add_comparison(
        session,
        configuration_id=config.configuration_id,
        query_hash=query_md5,
        subject_hash=subject_md5,
        identity=identity,
        aln_length=aln_length,
        sim_errors=sim_errors,
    )
    session.commit()
    return new


@app.command()
def log_comparison(  # noqa: PLR0913
    database: Annotated[str, typer.Option(help="Path to pyANI-plus SQLite3 database")],
    # These are for the comparison table
    query_fasta: Annotated[Path, typer.Option(help="Path to query FASTA file")],
    subject_fasta: Annotated[Path, typer.Option(help="Path to subject FASTA file")],
    identity: Annotated[
        float, typer.Option(help="Percent identity (float from 0 to 1)")
    ],
    aln_length: Annotated[int, typer.Option(help="Alignment length")],
    # These are all for the configuration table:
    method: Annotated[str, typer.Option(help="Comparison method")],
    program: Annotated[str, typer.Option(help="Comparison program name")],
    version: Annotated[str, typer.Option(help="Comparison program version")],
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
) -> int:
    """Log a single pyANI-plus pairwise comparison to the database."""
    print(f"Logging to {database}")  # noqa: T201
    db = Path(database)
    if not db.is_file():
        sys.exit(f"ERROR - Database {db} does not exist")
    session = db_orm.connect_to_db(db)

    config = db_orm.add_configuration(
        session, method, program, version, fragsize, maxmatch, kmersize, minmatch
    )
    if config.configuration_id is None:
        session.commit()
        print(f"Adding novel configuration {config.configuration_id} to database")  # noqa: T201
    else:
        print(f"Using pre-existing configuration {config.configuration_id} in database")  # noqa: T201

    new = _log_comparison(
        session, config, query_fasta, subject_fasta, identity, aln_length
    )
    if new:
        print(f"{query_fasta} vs {subject_fasta} logged")  # noqa: T201
    else:
        print(f"{query_fasta} vs {subject_fasta} existed")  # noqa: T201
    return 0


def parse_fastani_file(filename: Path) -> tuple[Path, Path, float, int, int]:
    """Parse a single-line fastANI output file extracting key fields as a tuple.

    Return (ref genome, query genome, ANI estimate, orthologous matches,
    sequence fragments) tuple.

    :param filename: Path, path to the input file

    Extracts the ANI estimate, the number of orthologous matches, and the
    number of sequence fragments considered from the fastANI output file.

    We assume that all fastANI comparisons are pairwise: one query and
    one reference file. The fastANI file should contain a single line.

    fsatANI *can* produce multi-line output, if a list of query/reference
    files is given to it.
    """
    with filename.open() as handle:
        line = handle.readline().strip().split()
    if not line:  # No file content; either run failed or no detectable similarity
        msg = f"Input file {filename} is empty"
        raise ValueError(msg)
    return (
        Path(line[0]),
        Path(line[1]),
        0.01 * float(line[2]),
        int(line[3]),
        int(line[4]),
    )


@app.command()
def log_fastani(  # noqa: PLR0913
    database: Annotated[str, typer.Option(help="Path to pyANI-plus SQLite3 database")],
    # These are for the comparison table
    query_fasta: Annotated[Path, typer.Option(help="Path to query FASTA file")],
    subject_fasta: Annotated[Path, typer.Option(help="Path to subject FASTA file")],
    fastani: Annotated[Path, typer.Option(help="Path to fastANI output file")],
    # These are all for the configuration table:
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
) -> int:
    """Log a single pyANI-plus fastANI pairwise comparison to the database."""
    # Assuming this will match as expect this script to be called right
    # after the computation has finished (on the same machine)
    fastani_tool = tools.get_fastani()

    used_query, used_subject, identity, aln_length, fragments = parse_fastani_file(
        fastani
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

    print(f"Logging to {database}")  # noqa: T201
    db = Path(database)
    if not db.is_file():
        sys.exit(f"ERROR - Database {db} does not exist")
    session = db_orm.connect_to_db(db)

    config = db_orm.add_configuration(
        session,
        method="fastANI",
        program=fastani_tool.exe_path.stem,
        version=fastani_tool.version,
        fragsize=fragsize,
        maxmatch=maxmatch,
        kmersize=kmersize,
        minmatch=minmatch,
    )
    if config.configuration_id is None:
        session.commit()
        print(f"Adding novel configuration {config.configuration_id} to database")  # noqa: T201
    else:
        print(f"Using pre-existing configuration {config.configuration_id} in database")  # noqa: T201

    sim_errors = fragments - aln_length  # here aln_length was fastANI's total length
    # Need to lookup query length for: query_cover = float(tot_length) / org_lengths[qname]

    new = _log_comparison(
        session,
        config,
        query_fasta=query_fasta,
        subject_fasta=subject_fasta,
        identity=identity,
        aln_length=aln_length,
        sim_errors=sim_errors,
    )
    if new:
        print(f"{query_fasta} vs {subject_fasta} logged")  # noqa: T201
    else:
        print(f"{query_fasta} vs {subject_fasta} existed")  # noqa: T201
    return 0


if __name__ == "__main__":
    sys.exit(app())
