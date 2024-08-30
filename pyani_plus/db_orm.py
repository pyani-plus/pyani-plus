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

from sqlalchemy import (
    Boolean,
    Column,
    Float,
    ForeignKey,
    Integer,
    String,
    UniqueConstraint,
    create_engine,
)
from sqlalchemy.orm import DeclarativeBase, Session, relationship, sessionmaker


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

    query_comparisons = relationship(
        "Comparison",
        back_populates="query",
        primaryjoin="Genome.genome_hash == Comparison.query_hash",
    )
    subject_comparisons = relationship(
        "Comparison",
        back_populates="subject",
        primaryjoin="Genome.genome_hash == Comparison.subject_hash",
    )

    def __repr__(self) -> str:
        """Return string representation of Genome table object."""
        return (
            f"Genome(genome_hash={self.genome_hash!r}, path={self.path!r},"
            f" length={self.length}, description={self.description!r})"
        )


class Configuration(Base):
    """Describes the configuration of a pyANI-plus run.

    Each run is a set of all-vs-all comparisons between a set of FASTA
    files (genome table rows), for a specific ANI algorithm or method,
    recorded as a set of comparison table rows.

    Each Run entry represents a specific instance of an analysis run which
    records details like the exact command line, and when it was run. It will
    point to a single Configuration entry (this class) which records the tool
    name, version, and parameters like k-mer size (some parameters are tool
    specific).

    Comparison entries are cached for unique combinations of the query and
    subject genomes (via the hash of their file contents, not via filenames),
    AND the run configuration.
    """

    __tablename__ = "configuration"

    __table_args__ = (
        UniqueConstraint(
            "method",
            "program",
            "version",
            "fragsize",
            "maxmatch",
            "kmersize",
            "minmatch",
        ),
    )

    configuration_id = Column(Integer, primary_key=True)
    comparisons = relationship("Comparison", back_populates="configuration")

    # This was part of the Run table in pyANI v0.2
    method = Column(String)
    # These were all part of the Comparison table in pyANI v0.2, which had
    # them as part of the uniqueness constraint.
    program = Column(String)
    version = Column(String)
    fragsize = Column(Integer)  # in fastANI this is fragLength
    maxmatch = Column(Boolean)  # in fastANi this is Null
    kmersize = Column(Integer)
    minmatch = Column(Float)

    def __repr__(self) -> str:
        """Return string representation of Genome table object."""
        return (
            f"Configuration(configuration_id={self.configuration_id},"
            f" program={self.program!r}, version={self.version!r},"
            f" fragsize={self.fragsize}, maxmatch={self.maxmatch},"
            f" kmersize={self.kmersize}, minmatch={self.maxmatch})"
        )


class Comparison(Base):
    """Describes a single pairwise comparison between two genomes (query and subject)."""

    __tablename__ = "comparisons"
    __table_args__ = (
        UniqueConstraint(
            "query_hash",
            "subject_hash",
            "configuration_id",
        ),
    )

    comparison_id = Column(Integer, primary_key=True)

    # See https://docs.sqlalchemy.org/en/20/orm/basic_relationships.html
    query_hash = Column(String, ForeignKey("genomes.genome_hash"), nullable=False)
    query = relationship(Genome, foreign_keys=[query_hash])

    subject_hash = Column(String, ForeignKey("genomes.genome_hash"), nullable=False)
    subject = relationship(Genome, foreign_keys=[subject_hash])

    configuration_id = Column(
        Integer, ForeignKey("configuration.configuration_id"), nullable=False
    )
    configuration = relationship(Configuration, foreign_keys=[configuration_id])

    # The results of the comparison
    identity = Column(Float)
    aln_length = Column(Integer)  # in fastANI this is matchedfrags * fragLength
    sim_errs = Column(Integer)  # in fastANI this is allfrags - matchedfrags
    cov_query = Column(Float)  # in fastANI this is matchedfrags/allfrags
    cov_subject = Column(Float)  # in fastANI this is Null

    def __str__(self) -> str:
        """Return string summarising the Comparison table row."""
        return (
            f"Query: {self.query_hash}, Subject: {self.subject_hash}, "
            f"%%ID={self.identity}, ({self.configuration.program} {self.configuration.version}), "
            f"FragSize: {self.configuration.fragsize}, MaxMatch: {self.configuration.maxmatch}, "
            f"KmerSize: {self.configuration.kmersize}, MinMatch: {self.configuration.minmatch}"
        )

    def __repr__(self) -> str:
        """Return string representation of Comparison table object."""
        return (
            f"Comparison(comparison_id={self.comparison_id!r}, "
            f"query_hash={self.query_hash!r}, "
            f"subject_hash={self.subject_hash!r}, "
            f"configuration_id={self.configuration_id!r}, "
            f"identity={self.identity}, "
            f"aln_length={self.aln_length}, "
            f"sim_errs={self.sim_errs}, "
            f"cov_query={self.cov_query}, "
            f"cov_subject={self.cov_subject})"
        )


def connect_to_db(dbpath: Path, *, echo: bool = False) -> Session:
    """Create/connect to existing DB, and return session bound to it.

    >>> session = connect_to_db("/tmp/example.sqlite", echo=True)
    """
    engine = create_engine(url=f"sqlite:///{dbpath}", echo=echo)
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine)()
