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
from Bio.SeqIO.FastaIO import SimpleFastaParser
from rich.progress import track

from pyani_plus import db_orm
from pyani_plus.utils import file_md5sum

app = typer.Typer()


@app.command()
def log_genome(
    database: Annotated[str, typer.Option(help="Path to pyANI-plus SQLite3 database")],
    fasta: Annotated[
        list[str] | None, typer.Argument(help="Path to FASTA file(s)")
    ] = None,
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
            if session.query(db_orm.Genome).where(db_orm.Genome.genome_hash == md5):
                file_skip += 1
                continue
            length = 0
            description = None
            with Path(filename).open() as handle:
                for title, seq in SimpleFastaParser(handle):
                    length += len(seq)
                    if description is None:
                        description = title  # Just use first entry
            genome = db_orm.Genome(
                genome_hash=md5, path=filename, length=length, description=description
            )
            session.add(genome)
    session.commit()
    session.close()
    print(  # noqa: T201
        f"Processed {file_total} FASTA files, skipped {file_skip}, recorded {file_total-file_skip}"
    )

    return 0


@app.command()
def log_comparison(
    query_fasta: Annotated[str, typer.Argument(help="Path to query FASTA file")],
    subject_fasta: Annotated[str, typer.Argument(help="Path to subject FASTA file")],
    database: Annotated[str, typer.Option(help="Path to pyANI-plus SQLite3 database")],
) -> int:
    """Log a single pyANI-plus pairwise comparison to the database."""
    print(f"Logging to {database}")  # noqa: T201
    print(f"{query_fasta} vs {subject_fasta}")  # noqa: T201
    return 0


if __name__ == "__main__":
    sys.exit(app())
