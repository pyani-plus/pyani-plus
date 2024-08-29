# Copyright (c) 2024-present University of Strathclyde
# Author: Peter Cock
#
# Contact:
# peter.cock@strath.ac.uk
#
# Peter Cock,
# Strathclyde Institute of Pharmaceutical and Biomedical Sciences
# The University of Strathclyde
# 161 Cathedral Street
# Glasgow
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2024-present University of Strathclyde
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
"""Object Relational Mapping for a pyANI-plus SQLite3 database.

Using SQLalchemy, the Python classes defined here give us a
database schema and the code to import/export the data as
Python objects.
"""

from pathlib import Path

from sqlalchemy import Column, Integer, String, create_engine
from sqlalchemy.orm import DeclarativeBase, Session, sessionmaker


class Base(DeclarativeBase):
    """Base class for SQLAlchemy ORM declarations.

    See the SQLAlchemy 2.0 documentation. This is expected to be
    compatible with type checkers like mypy.
    """


class Genome(Base):
    """Describes an input genome (FASTA file) for a pyANI-plus run.

    We identify genome FASTA files by the MD5 checksum of the file contents.
    This allows us to treat copies of the same FASTA file the same, even
    if renamed or moved. This does not extend to any edits to the file,
    even if the sequence is unchanged.

    - genome_hash
        primary key, MD5 hash of input genome file (in ``path``)
    - path
        path to FASTA genome file (as first seen by pyANI)
    - length
        length of genome (total bases)
    - description
        genome description (defaults to the description of the first sequence)
    """

    __tablename__ = "genomes"

    genome_hash = Column(String, primary_key=True)
    path = Column(String)
    length = Column(Integer)  # total length of all the sequences
    description = Column(String)

    def __repr__(self) -> str:
        """Return string representation of Genome table object."""
        return (
            f"Genome(genome_hash={self.genome_hash!r}, path={self.path!r},"
            f" length={self.length}, description={self.description!r})"
        )


def connect_to_db(dbpath: Path, *, echo: bool = False) -> Session:
    """Create/connect to existing DB, and return session bound to it.

    >>> session = connect_to_db("/tmp/example.sqlite", echo=True)
    """
    engine = create_engine(url=f"sqlite:///{dbpath}", echo=echo)
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine)()
