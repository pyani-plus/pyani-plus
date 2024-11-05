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
import os
import platform
from collections.abc import Callable
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from sqlalchemy import (
    Dialect,
    ForeignKey,
    String,
    UniqueConstraint,
    create_engine,
    insert,
)
from sqlalchemy.exc import NoResultFound
from sqlalchemy.orm import (
    DeclarativeBase,
    Mapped,
    Session,
    aliased,
    mapped_column,
    object_session,
    registry,
    relationship,
    sessionmaker,
)
from sqlalchemy.types import TypeDecorator


class PathLike(TypeDecorator):
    """Allows mapping an `os.PathLike` object to an `sqlalchemy.String`.

    This can be used with `pathlib.Path` and other subclasses, as long
    as they can be turned into strings.
    """

    impl = String

    def __init__(self, factory: Callable[[str], os.PathLike]) -> None:
        """Create a sqlalchemy mapping between Python paths and SQL strings."""
        super().__init__()
        self.factory = factory

    def process_bind_param(
        self,
        value: os.PathLike | None,
        dialect: Dialect,  # noqa: ARG002
    ) -> str | None:
        """Convert an `os.PathLike` value to a string for the database."""
        # What if the value is bytes not a string?
        return os.fspath(value) if value else None

    def process_result_value(
        self,
        value: str | None,
        dialect: Dialect,  # noqa: ARG002
    ) -> os.PathLike | None:
        """Restore a string from the database to an `os.Pathlike`."""
        return None if value is None else self.factory(value)


class Base(DeclarativeBase):
    """Base class for SQLAlchemy ORM declarations.

    See the SQLAlchemy 2.0 documentation. This is expected to be
    compatible with type checkers like mypy.
    """

    registry = registry(type_annotation_map={Path: PathLike(Path)})


class RunGenomeAssociation(Base):
    """Linker table between genomes and runs tables."""

    __tablename__ = "runs_genomes"
    genome_hash: Mapped[str] = mapped_column(
        ForeignKey("genomes.genome_hash"), primary_key=True
    )
    run_id: Mapped[int] = mapped_column(ForeignKey("runs.run_id"), primary_key=True)
    # Just filename, no containing directory - that's in the associated run table entry
    fasta_filename: Mapped[str] = mapped_column()

    run = relationship("Run")
    genome = relationship("Genome")


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
    path: Mapped[Path] = mapped_column()
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
    runs = relationship("Run", secondary="runs_genomes", lazy="dynamic", viewonly=True)

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

    __tablename__ = "configurations"

    __table_args__ = (
        UniqueConstraint(
            "method",
            "program",
            "version",
            "fragsize",
            "mode",
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

    # This was part of the Run table in pyANI v0.2
    method: Mapped[str] = mapped_column()
    # These were all part of the Comparison table in pyANI v0.2, which had
    # them as part of the uniqueness constraint.
    program: Mapped[str] = mapped_column()
    version: Mapped[str] = mapped_column()
    fragsize: Mapped[int | None] = mapped_column()  # in fastANI this is fragLength
    mode: Mapped[str | None] = mapped_column()  # this is "mum" or "maxmatch" in ANIm
    kmersize: Mapped[int | None] = mapped_column()
    minmatch: Mapped[float | None] = mapped_column()

    def __repr__(self) -> str:
        """Return string representation of Configuration table object."""
        return (
            f"Configuration(configuration_id={self.configuration_id},"
            f" program={self.program!r}, version={self.version!r},"
            f" fragsize={self.fragsize}, mode={self.mode!r},"
            f" kmersize={self.kmersize}, minmatch={self.minmatch})"
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
        ForeignKey("configurations.configuration_id"), nullable=False
    )
    configuration: Mapped[Configuration] = relationship(
        Configuration, foreign_keys=[configuration_id]
    )

    # The results of the comparison
    identity: Mapped[float] = mapped_column()
    # in fastANI this is matchedfrags * fragLength:
    aln_length: Mapped[int] = mapped_column()
    # in fastANI this is allfrags - matchedfrags
    sim_errors: Mapped[int | None] = mapped_column()
    # in fastANI this is matchedfrags/allfrags
    cov_query: Mapped[float | None] = mapped_column()
    # in fastANI this is Null
    cov_subject: Mapped[float | None] = mapped_column()

    # Logging about where the comparison was run (which is deliberately not
    # part of the configuration class). Unlike underlying tool versions, we
    # don't expect the OS and architecture to matter - but they could do,
    # most likely though minor floating point differences. By logging this
    # here, we can more easily compare the same run on different platforms
    # (each run with its own local database):
    #
    # e.g. "Darwin" from `uname -s` or platform.uname().system
    uname_system: Mapped[str] = mapped_column()
    # e.g. "21.6.0" from `uname -r` or platform.uname().release
    uname_release: Mapped[str] = mapped_column()
    # e.g. "arm64" from `uname -m` or platform.uname().machine
    uname_machine: Mapped[str] = mapped_column()

    def __str__(self) -> str:
        """Return string summarising the Comparison table row."""
        return (
            f"Query: {self.query_hash}, Subject: {self.subject_hash}, "
            f"%ID={self.identity}, ({self.configuration.program} {self.configuration.version}), "
            f"FragSize: {self.configuration.fragsize}, Mode: {self.configuration.mode}, "
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
            f"sim_errors={self.sim_errors}, "
            f"cov_query={self.cov_query}, "
            f"cov_subject={self.cov_subject}, "
            f"uname_system={self.uname_system!r}, "
            f"uname_release={self.uname_release!r}, "
            f"uname_machine={self.uname_machine!r})"
        )


class Run(Base):
    """Describes a single pyANI-plus run.

    Each run is a set of all-vs-all comparisons between a set of FASTA
    files (genome table rows), for a specific ANI algorithm or method
    (represented as a row in the configuration table), recorded as a set
    of comparison table rows.

    Thus one Run is linked to one Configuration row, many Genome rows,
    and many Comparison rows, where each Comparison is linked to two
    Genome rows (query and reference).

    If any of these these comparisons have been computed and recorded by
    a previous run using the same configuration and an overlapping set
    of genomes, the new run will point to those pre-existing comparisons.

    In order to generate reports and plots quickly, the Run entry also
    caches dataframes of the information in the linked Comparison entries
    (as JSON-encoded Pandas dataframes).
    """

    __tablename__ = "runs"

    run_id: Mapped[int] = mapped_column(primary_key=True)

    configuration_id: Mapped[int] = mapped_column(
        ForeignKey("configurations.configuration_id"), nullable=False
    )
    configuration: Mapped[Configuration] = relationship(
        Configuration, foreign_keys=[configuration_id]
    )

    cmdline: Mapped[str] = mapped_column()
    fasta_directory: Mapped[Path] = mapped_column()
    date: Mapped[datetime.datetime] = mapped_column()
    status: Mapped[str] = mapped_column()
    name: Mapped[str] = mapped_column()
    df_identity: Mapped[str | None] = mapped_column()  # JSON-encoded Pandas dataframe
    df_cov_query: Mapped[str | None] = mapped_column()  # JSON-encoded Pandas dataframe
    df_aln_length: Mapped[str | None] = mapped_column()  # JSON-encoded Pandas dataframe
    df_sim_errors: Mapped[str | None] = mapped_column()  # JSON-encoded Pandas dataframe
    df_hadamard: Mapped[str | None] = mapped_column()  # JSON-encoded Pandas dataframe

    fasta_hashes: Mapped[list[RunGenomeAssociation]] = relationship(
        RunGenomeAssociation, viewonly=True, lazy="dynamic"
    )
    genomes: Mapped[list[Genome]] = relationship(
        Genome, secondary="runs_genomes", viewonly=True, lazy="dynamic"
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
            .where(Comparison.configuration_id == self.configuration_id)
            .join(run_query, Comparison.query_hash == run_query.genome_hash)
            .join(run_subjt, Comparison.subject_hash == run_subjt.genome_hash)
            .where(run_query.run_id == self.run_id)
            .where(run_subjt.run_id == self.run_id)
        )

    def cache_comparisons(self) -> pd.DataFrame:
        """Collect and cache all-vs-all score matrices for the run.

        If the run has N genomes, should have an N by N matrix of all the pairwise
        comparisons, coming from N^2 entries in the comparisons table. These
        will be fetched, converted into pandas dataframes, stored in the attributes
        .df_identity, .df_coverage, .df_aln_length, .df_sim_errors, and .df_hadamard
        as JSON string using the compact "split" orientation.

        The caller must commit the updated Run object to the database explicitly!
        """
        hashes = sorted(association.genome_hash for association in self.fasta_hashes)
        size = len(hashes)
        identity = np.full([size, size], np.nan, float)
        cov_query = np.full([size, size], np.nan, float)
        aln_length = np.full([size, size], np.nan, float)
        sim_errors = np.full([size, size], np.nan, float)
        for comp in self.comparisons():
            row = hashes.index(comp.query_hash)
            col = hashes.index(comp.subject_hash)
            identity[row, col] = comp.identity
            cov_query[row, col] = comp.cov_query
            aln_length[row, col] = comp.aln_length
            sim_errors[row, col] = comp.sim_errors
        # Hadamard matrix is (element wise) identity * coverage
        self.df_hadamard = pd.DataFrame(
            data=identity * cov_query, index=hashes, columns=hashes, dtype=float
        ).to_json(orient="split")
        self.df_identity = pd.DataFrame(
            data=identity, index=hashes, columns=hashes, dtype=float
        ).to_json(orient="split")
        del identity
        self.df_cov_query = pd.DataFrame(
            data=cov_query, index=hashes, columns=hashes, dtype=float
        ).to_json(orient="split")
        del cov_query
        self.df_aln_length = pd.DataFrame(
            data=aln_length, index=hashes, columns=hashes, dtype=float
        ).to_json(orient="split")
        del aln_length
        self.df_sim_errors = pd.DataFrame(
            data=sim_errors, index=hashes, columns=hashes, dtype=float
        ).to_json(orient="split")
        del sim_errors

    @property
    def identities(self) -> pd.DataFrame | None:
        """All-vs-all percentage identity matrix for the run from the cached JSON.

        If cached, returns an N by N float matrix of percentage identities for the N genomes
        in the run as pandas dataframe, where the index (rows) and columns are the N genome
        hashes (sorted alphabetically). If not cached, returns None.

        This is normally available immediately from the Run entry in the database where
        the dataframe is cached as a JSON string. This would be computed at the end of
        the run once all N^2 comparisons have finished, see the cache_comparisons method.
        """
        if not self.df_identity:
            return None
        return pd.read_json(StringIO(self.df_identity), orient="split", dtype=float)

    @property
    def cov_query(self) -> pd.DataFrame | None:
        """All-vs-all query-coverage matrix for the run from the cached JSON.

        If cached, returns an N by N float matrix of percentage identities for the N genomes
        in the run as a pandas dataframe, where the index (rows) and columns are the N
        genome hashes (sorted alphabetically). If not cached, returns None.

        This is normally available immediately from the Run entry in the database where
        the dataframe is cached as a JSON string. This would be computed at the end of
        the run once all N^2 comparisons have finished, see the cache_comparisons method.
        """
        if not self.df_cov_query:
            return None
        return pd.read_json(StringIO(self.df_cov_query), orient="split", dtype=float)

    @property
    def aln_length(self) -> pd.DataFrame | None:
        """All-vs-all alignment length matrix for the run from the cached JSON.

        If cached, returns an N by N integer matrix of alignment lengths for the N genomes
        in the run as a pandas dataframe, where the index (rows) and columns are the N
        genome hashes (sorted alphabetically). If not cached, returns None.

        This is normally available immediately from the Run entry in the database where
        the dataframe is cached as a JSON string. This would be computed at the end of
        the run once all N^2 comparisons have finished, see the cache_comparisons method.
        """
        if not self.df_aln_length:
            return None
        return pd.read_json(StringIO(self.df_aln_length), orient="split")

    @property
    def sim_errors(self) -> pd.DataFrame | None:
        """All-vs-all similarity errors matrix for the run from the cached JSON.

        If cached, returns an N by N integer matrix of similarity errors for the N genomes
        in the run as a pandas dataframe, where the index (rows) and columns are the N
        genome hashes (sorted alphabetically). If not cached, returns None.

        This is normally available immediately from the Run entry in the database where
        the dataframe is cached as a JSON string. This would be computed at the end of
        the run once all N^2 comparisons have finished, see the cache_comparisons method.
        """
        if not self.df_sim_errors:
            return None
        return pd.read_json(StringIO(self.df_sim_errors), orient="split")

    @property
    def hadamard(self) -> pd.DataFrame | None:
        """All-vs-all Hadamard matrix (identity times coverage) for the run from the cached JSON.

        If cached, returns an N by N Hadamard matrix for the N genomes in the run as a
        pandas dataframe, where the index (rows) and columns are the N genome hashes
        (sorted alphabetically). If not cached, returns None.

        This is normally available immediately from the Run entry in the database where
        the dataframe is cached as a JSON string. This would be computed at the end of
        the run once all N^2 comparisons have finished, see the cache_comparisons method.
        """
        # Might be worth profiling the performance of caching this, versus
        # computing it from the cached identity and coverage data-frames?
        if not self.df_hadamard:
            return None
        return pd.read_json(StringIO(self.df_hadamard), orient="split", dtype=float)

    def __repr__(self) -> str:
        """Return abridged string representation of Run table object."""
        return (
            f"Run(run_id={self.run_id}, configuration_id={self.configuration_id},"
            f" cmdline={self.cmdline!r}, date={self.date!r},"
            f" status={self.status!r}, name={self.name!r}, ...)"
        )


def connect_to_db(dbpath: Path | str, *, echo: bool = False) -> Session:
    """Create/connect to existing DB, and return session bound to it.

    >>> session = connect_to_db("/tmp/pyani-plus-example.sqlite", echo=True)
    20...

    Will accept the special SQLite3 value of ":memory:" for an in-memory
    database:

    >>> session = connect_to_db(":memory:")
    """
    # Note with echo=True, the output starts yyyy-mm-dd and sadly
    # using just ... is interpreted as a continuation of the >>>
    # prompt rather than saying any output is fine with ELLIPSIS mode.

    # Default timeout is 5s
    engine = create_engine(
        url=f"sqlite:///{dbpath!s}", echo=echo, connect_args={"timeout": 10}
    )
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine)()


def db_configuration(  # noqa: PLR0913
    session: Session,
    method: str,
    program: str,
    version: str,
    fragsize: int | None = None,
    mode: str | None = None,
    kmersize: int | None = None,
    minmatch: float | None = None,
    *,
    create: bool = False,
) -> Configuration:
    """Return a configuration table entry, or add and return it if not already there.

    By default if the entry is not there already, you get a NoResultFound exception:

    >>> session = connect_to_db(":memory:")
    >>> conf = db_configuration(
    ...     session,
    ...     method="guessing",
    ...     program="guestimate",
    ...     version="v0.1.2beta3",
    ...     fragsize=1000,
    ...     kmersize=31,
    ... )
    Traceback (most recent call last):
    ...
    sqlalchemy.exc.NoResultFound: Requested configuration not already in DB

    If the entry is not there already, and you want to add it, you must use create=True:

    >>> session = connect_to_db(":memory:")
    >>> conf = db_configuration(
    ...     session,
    ...     method="guessing",
    ...     program="guestimate",
    ...     version="v0.1.2beta3",
    ...     fragsize=1000,
    ...     kmersize=31,
    ...     create=True,
    ... )
    >>> conf.configuration_id
    1

    Note this calls session.commit() explicitly to try to reduce locking contention.
    """
    config = (
        session.query(Configuration)
        .where(Configuration.method == method)
        .where(Configuration.program == program)
        .where(Configuration.version == version)
        .where(Configuration.fragsize == fragsize)
        .where(Configuration.mode == mode)
        .where(Configuration.kmersize == kmersize)
        .where(Configuration.minmatch == minmatch)
        .one_or_none()
    )
    if config is None:
        if not create:
            msg = "Requested configuration not already in DB"
            raise NoResultFound(msg)
        config = Configuration(
            method=method,
            program=program,
            version=version,
            fragsize=fragsize,
            mode=mode,
            kmersize=kmersize,
            minmatch=minmatch,
        )
        session.add(config)
        session.commit()
    return config


def db_genome(
    session: Session, fasta_filename: Path, md5: str, *, create: bool = False
) -> Genome:
    """Return a genome table entry, or add and return it if not already there.

    Assumes and trusts the MD5 checksum given matches.

    Returns the matching genome object, or the new one added if create=True:

    >>> session = connect_to_db(":memory:")
    >>> from pyani_plus.utils import file_md5sum
    >>> fasta = Path("tests/fixtures/viral_example/OP073605.fasta")
    >>> genome = db_genome(session, fasta, file_md5sum(fasta), create=True)
    >>> genome.genome_hash
    '5584c7029328dc48d33f95f0a78f7e57'

    Note this calls session.commit() explicitly to try to reduce locking contention,
    and makes some efforts to cope gracefully with competing processes also trying
    to record the same genome at the same time.

    If the genome is not already there, then by default this raises an exception:

    >>> session = connect_to_db(":memory:")
    >>> genome = db_genome(session, fasta, file_md5sum(fasta))
    Traceback (most recent call last):
    ...
    sqlalchemy.exc.NoResultFound: Requested genome not already in DB
    """
    old_genome = session.query(Genome).where(Genome.genome_hash == md5).one_or_none()
    if old_genome is not None:
        return old_genome

    if not create:
        msg = "Requested genome not already in DB"
        raise NoResultFound(msg)

    length = 0
    description = None
    with fasta_filename.open() as handle:
        for title, seq in SimpleFastaParser(handle):
            length += len(seq)
            if description is None:
                description = title  # Just use first entry
    genome = Genome(
        genome_hash=md5,
        path=fasta_filename,
        length=length,
        description=description,
    )

    # Recheck database as all that file access is slow and another thread
    # may have added it in the meantime
    old_genome = session.query(Genome).where(Genome.genome_hash == md5).one_or_none()
    if old_genome is not None:
        return old_genome

    session.add(genome)
    session.commit()
    return genome


def add_run(  # noqa: PLR0913
    session: Session,
    configuration: Configuration,
    cmdline: str,
    fasta_directory: Path,
    status: str,
    name: str,
    date: datetime.datetime | None = None,
    fasta_to_hash: dict[Path, str] | None = None,
) -> Run:
    """Add and return a new run table entry.

    Will make a near-duplicate if there is a match there already!

    """
    run = Run(
        configuration_id=configuration.configuration_id,
        cmdline=cmdline,
        fasta_directory=fasta_directory,
        status=status,
        name=name,
        date=date if date else datetime.datetime.now(tz=datetime.UTC),
    )
    session.add(run)
    if fasta_to_hash:
        # Using the low-level but faster insert we need the run_id, so must commit now:
        session.commit()
        run_id = run.run_id
        session.execute(
            insert(RunGenomeAssociation),
            [
                # Note /mnt/data/example.fasta becomes just example.fasta
                # with the path /mnt/data/ recorded in the run itself
                # (since this is by design shared for all FASTA in a run)
                {"run_id": run_id, "fasta_filename": filename.name, "genome_hash": md5}
                for filename, md5 in fasta_to_hash.items()
            ],
        )
    session.commit()
    return run


def db_comparison(  # noqa: PLR0913
    session: Session,
    configuration_id: int,
    query_hash: str,
    subject_hash: str,
    identity: float,
    aln_length: int,
    sim_errors: int | None = None,
    cov_query: float | None = None,
    cov_subject: float | None = None,
    uname: platform.uname_result | None = None,
) -> Comparison:
    """Return a comparison table entry, or add and return it if not already there.

    This assumes the configuration and both the query and subject are already in
    the linked tables. If not, addition will fail with an integrity error:

    >>> session = connect_to_db(":memory:")
    >>> comp = db_comparison(session, 1, "abcd", "cdef", 0.99, 12345)
    >>> comp.identity
    0.99

    By default the uname values for the current platform are used (i.e. this assumes
    the computed values were computed locally on the same machine).

    This will NOT alter a pre-existing entry even if there are differences (e.g.
    expected like operating system, or unexpected like percentage identity).

    Note this calls session.commit() explicitly to try to reduce locking contention.
    """
    comp = (
        session.query(Comparison)
        .where(Comparison.configuration_id == configuration_id)
        .where(Comparison.query_hash == query_hash)
        .where(Comparison.subject_hash == subject_hash)
        .one_or_none()
    )
    if comp is not None:
        # We could sanity check the entry there matches... that might reveal
        # a cross-platform difference, or a change to historical behaviour?
        return comp

    if uname is None:
        # This function caches the return value, so repeat calls are fast:
        uname = platform.uname()

    comp = Comparison(
        configuration_id=configuration_id,
        query_hash=query_hash,
        subject_hash=subject_hash,
        identity=identity,
        aln_length=aln_length,
        sim_errors=sim_errors,
        cov_query=cov_query,
        cov_subject=cov_subject,
        uname_system=uname.system,
        uname_release=uname.release,
        uname_machine=uname.machine,
    )
    session.add(comp)
    session.commit()
    return comp
