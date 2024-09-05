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
"""Object Relational Mapping for a pyANI-plus SQLite3 database.

Using SQLalchemy, the Python classes defined here give us a
database schema and the code to import/export the data as
Python objects.
"""

import datetime
from pathlib import Path

from sqlalchemy import (
    ForeignKey,
    UniqueConstraint,
    create_engine,
)
from sqlalchemy.orm import (
    DeclarativeBase,
    Mapped,
    Session,
    aliased,
    mapped_column,
    object_session,
    relationship,
    sessionmaker,
)


class Base(DeclarativeBase):
    """Base class for SQLAlchemy ORM declarations.

    See the SQLAlchemy 2.0 documentation. This is expected to be
    compatible with type checkers like mypy.
    """


class RunGenomeAssociation(Base):
    """Linker table between genomes and runs tables."""

    __tablename__ = "runs_genomes"
    genome_hash: Mapped[str] = mapped_column(
        ForeignKey("genomes.genome_hash"), primary_key=True
    )
    run_id: Mapped[int] = mapped_column(ForeignKey("runs.run_id"), primary_key=True)


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

    genome_hash: Mapped[str] = mapped_column(primary_key=True)
    path: Mapped[str] = mapped_column()
    length: Mapped[int] = mapped_column()  # total length of all the sequences
    description: Mapped[str] = mapped_column()

    query_comparisons: Mapped[list["Comparison"]] = relationship(
        "Comparison",
        back_populates="query",
        primaryjoin="Genome.genome_hash == Comparison.query_hash",
    )
    subject_comparisons: Mapped[list["Comparison"]] = relationship(
        "Comparison",
        back_populates="subject",
        primaryjoin="Genome.genome_hash == Comparison.subject_hash",
    )
    runs = relationship(
        "Run", secondary="runs_genomes", back_populates="genomes", lazy="dynamic"
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

    configuration_id: Mapped[int] = mapped_column(primary_key=True)
    comparisons: Mapped[list["Comparison"]] = relationship(
        "Comparison",
        back_populates="configuration",
        primaryjoin="Configuration.configuration_id == Comparison.configuration_id",
        lazy="dynamic",
    )
    runs: Mapped[list["Run"]] = relationship(
        "Run",
        back_populates="configuration",
        primaryjoin="Configuration.configuration_id == Run.configuration_id",
        lazy="dynamic",
    )

    # These are system properties which may affect the results
    # (most likely in fine details of floating point computations)
    machine: Mapped[str] = (
        mapped_column()
    )  # e.g. "arm64" from `uname -m` or platform.uname().machine
    system: Mapped[str] = (
        mapped_column()
    )  # e.g. "Darwin" from `uname -s` or platform.uname().system

    # This was part of the Run table in pyANI v0.2
    method: Mapped[str] = mapped_column()
    # These were all part of the Comparison table in pyANI v0.2, which had
    # them as part of the uniqueness constraint.
    program: Mapped[str] = mapped_column()
    version: Mapped[str] = mapped_column()
    fragsize: Mapped[int | None] = mapped_column()  # in fastANI this is fragLength
    maxmatch: Mapped[bool | None] = mapped_column()  # in fastANi this is Null
    kmersize: Mapped[int | None] = mapped_column()
    minmatch: Mapped[float | None] = mapped_column()

    def __repr__(self) -> str:
        """Return string representation of Genome table object."""
        return (
            f"Configuration(configuration_id={self.configuration_id},"
            f" machine={self.machine!r}, system={self.system!r},"
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

    comparison_id: Mapped[int] = mapped_column(primary_key=True)

    # See https://docs.sqlalchemy.org/en/20/orm/basic_relationships.html
    query_hash: Mapped[str] = mapped_column(
        ForeignKey("genomes.genome_hash"), nullable=False
    )
    query: Mapped[Genome] = relationship(Genome, foreign_keys=[query_hash])

    subject_hash: Mapped[str] = mapped_column(
        ForeignKey("genomes.genome_hash"), nullable=False
    )
    subject: Mapped[Genome] = relationship(Genome, foreign_keys=[subject_hash])

    configuration_id: Mapped[int] = mapped_column(
        ForeignKey("configuration.configuration_id"), nullable=False
    )
    configuration: Mapped[Configuration] = relationship(
        Configuration, foreign_keys=[configuration_id]
    )

    # The results of the comparison
    identity: Mapped[float] = mapped_column()
    # in fastANI this is matchedfrags * fragLength:
    aln_length: Mapped[int] = mapped_column()
    # in fastANI this is allfrags - matchedfrags
    sim_errs: Mapped[int | None] = mapped_column()
    # in fastANI this is matchedfrags/allfrags
    cov_query: Mapped[float | None] = mapped_column()
    # in fastANI this is Null
    cov_subject: Mapped[float | None] = mapped_column()

    def __str__(self) -> str:
        """Return string summarising the Comparison table row."""
        return (
            f"Query: {self.query_hash}, Subject: {self.subject_hash}, "
            f"%ID={self.identity}, ({self.configuration.program} {self.configuration.version}), "
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


class Run(Base):
    """Describes a single pyANI-plus run.

    Each run is a set of all-vs-all comparisons between a set of FASTA
    files (genome table rows), for a specific ANI algorithm or method
    (represented as a row in the configuration table),
    recorded as a set of comparison table rows.

    Thus one Run is linked to one Configuration row, many Genome rows,
    and many Comparison rows, where each Comparison is linked to two
    Genome rows (query and reference).

    In order to generate reports and plots quickly, the Run entry also
    caches dataframes of the information in the linked Comparison entries
    (as JSON-encoded Pandas dataframes).
    """

    __tablename__ = "runs"

    run_id: Mapped[int] = mapped_column(primary_key=True)

    configuration_id: Mapped[int] = mapped_column(
        ForeignKey("configuration.configuration_id"), nullable=False
    )
    configuration: Mapped[Configuration] = relationship(
        Configuration, foreign_keys=[configuration_id]
    )

    cmdline: Mapped[str] = mapped_column()
    date: Mapped[datetime.datetime] = mapped_column()
    status: Mapped[str] = mapped_column()
    name: Mapped[str] = mapped_column()
    df_identity: Mapped[str | None] = mapped_column()  # JSON-encoded Pandas dataframe
    df_coverage: Mapped[str | None] = mapped_column()  # JSON-encoded Pandas dataframe
    df_alnlength: Mapped[str | None] = mapped_column()  # JSON-encoded Pandas dataframe
    df_simerrors: Mapped[str | None] = mapped_column()  # JSON-encoded Pandas dataframe
    df_hadamard: Mapped[str | None] = mapped_column()  # JSON-encoded Pandas dataframe

    genomes: Mapped[list[Genome]] = relationship(
        Genome, secondary="runs_genomes", back_populates="runs", lazy="dynamic"
    )

    def comparisons(self) -> Mapped[list[Comparison]]:
        """Find all the comparison rows for this run.

        Runs a complex double join for the query and subject genomes
        of a comparison via the runs_genomes association table (using
        two aliases)::

            SELECT ... FROM comparisons
            JOIN runs_genomes AS run_query
            ON comparisons.query_hash = run_query.genome_hash
            JOIN runs_genomes AS run_subject
            ON comparisons.subject_hash = run_subject.genome_hash
            WHERE run_query.run_id = ? AND run_subject.run_id = ?

        This method is used internally as part of populating the cached
        identities, alignment lengths, etc for the run.
        """
        # See https://docs.sqlalchemy.org/en/20/orm/mapped_sql_expr.html
        # I couldn't work out how to do this with column_property etc.
        # It works as a simple property, but as an expensive search seems
        # clearer to define this as a class method instead.

        # This works using aliased(db_orm.RunGenomeAssociation) which
        # is a class-based definition of the linker table, but failed
        # when it was just aliased(rungenome) defined using
        # Table("runs_genomes", Base.metadata, ...) giving:
        # AttributeError: 'Alias' object has no attribute 'genome_hash'
        run_query = aliased(RunGenomeAssociation, name="run_query")
        run_subjt = aliased(RunGenomeAssociation, name="run_subject")

        return (
            object_session(self)
            .query(Comparison)
            .join(run_query, Comparison.query_hash == run_query.genome_hash)
            .join(run_subjt, Comparison.subject_hash == run_subjt.genome_hash)
            .where(run_query.run_id == self.run_id)
            .where(run_subjt.run_id == self.run_id)
        )

    def __repr__(self) -> str:
        """Return abridged string representation of Run table object."""
        return (
            f"Run(run_id={self.run_id}, configuration_id={self.configuration_id},"
            f" cmdline={self.cmdline!r}, date={self.date!r},"
            f" status={self.status!r}, name={self.name!r}, ...)"
        )


def connect_to_db(dbpath: Path, *, echo: bool = False) -> Session:
    """Create/connect to existing DB, and return session bound to it.

    >>> session = connect_to_db("/tmp/example.sqlite", echo=True)
    """
    engine = create_engine(url=f"sqlite:///{dbpath}", echo=echo)
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine)()
