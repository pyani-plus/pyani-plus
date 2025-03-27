# The MIT License
#
# Copyright (c) 2024-2025 University of Strathclyde
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
"""Tests for the pyani_plus/db_orm.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

import datetime
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sqlalchemy.exc import NoResultFound

from pyani_plus import db_orm, setup_logger
from pyani_plus.utils import file_md5sum, str_md5sum


def test_make_new_db(tmp_path: str) -> None:
    """Confirm can create a new empty database."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()
    logger = setup_logger(None)
    db_orm.connect_to_db(logger, tmp_db)  # discard the session, should close
    assert tmp_db.is_file()
    with tmp_db.open("rb") as handle:
        magic = handle.read(16)
        assert magic == b"SQLite format 3\0"
    tmp_db.unlink()


def test_make_and_populate_comparisons(tmp_path: str) -> None:
    """Populate new DB with config, genomes, and comparisons."""
    DUMMY_ALIGN_LEN = 400  # noqa: N806
    NAMES = ("Alpha", "Beta", "Gamma")  # noqa: N806
    tmp_db = Path(tmp_path) / "genomes.sqlite"
    assert not tmp_db.is_file()

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)

    uname = platform.uname()
    config = db_orm.Configuration(
        method="guessing",
        program="guestimate",
        version="v0.1.2beta3",
        fragsize=100,
        mode="RNG",
        kmersize=17,
        minmatch=0.3,
    )
    # Test the __repr__
    assert repr(config) == (
        "Configuration(configuration_id=None,"
        " program='guestimate', version='v0.1.2beta3',"
        " fragsize=100, mode='RNG', kmersize=17, minmatch=0.3, extra=None)"
    )
    session.add(config)
    assert config.configuration_id is None
    session.commit()
    assert config.configuration_id == 1

    hashes = {}
    for name in NAMES:
        seq = "ACGT" * int(DUMMY_ALIGN_LEN / 4) * len(name)
        fasta = f">{name}\n{seq}\n"
        md5 = str_md5sum(fasta)
        hashes[md5] = name
        genome = db_orm.Genome(
            genome_hash=md5,
            path=f"/mnt/shared/data/{name}.fasta",
            length=len(seq),
            description=f"Example {name}",
        )
        # Here testing the Genome class __repr__:
        assert repr(genome) == (
            f"Genome(genome_hash={md5!r}, path='/mnt/shared/data/{name}.fasta',"
            f" length={len(seq)}, description='Example {name}')"
        )
        session.add(genome)
    assert len(list(config.comparisons)) == 0
    assert session.query(db_orm.Genome).count() == len(NAMES)

    assert config.configuration_id == 1
    for query in hashes:
        for subject in hashes:
            comparison = db_orm.Comparison(
                # configuration=config,  <-- should this work?
                configuration_id=config.configuration_id,
                query_hash=query,
                subject_hash=subject,
                identity=0.96,
                aln_length=DUMMY_ALIGN_LEN,
                uname_machine=uname.machine,  # CPU arch
                uname_system=uname.system,  # Operating system
                uname_release=uname.release,
            )
            assert comparison.configuration_id == config.configuration_id
            assert repr(comparison) == (
                "Comparison(comparison_id=None,"
                f" query_hash={query!r}, subject_hash={subject!r},"
                f" configuration_id={config.configuration_id},"
                f" identity=0.96, aln_length={DUMMY_ALIGN_LEN}, sim_errors=None,"
                f" cov_query=None, cov_subject=None, uname_system={uname.system!r},"
                f" uname_release={uname.release!r}, uname_machine={uname.machine!r})"
            )
            # Can't test __str__ yet as .configuration not yet live?
            session.add(comparison)
    session.commit()
    assert session.query(db_orm.Comparison).count() == len(NAMES) ** 2
    assert len(list(config.comparisons)) == len(NAMES) ** 2
    for comparison in config.comparisons:
        assert comparison.aln_length == DUMMY_ALIGN_LEN
        assert str(comparison) == (
            f"Query: {comparison.query_hash}, Subject: {comparison.subject_hash},"
            " %ID=0.96, (guestimate v0.1.2beta3), "
            "FragSize: 100, Mode: RNG, KmerSize: 17, MinMatch: 0.3"
        )

        # Check the configuration object attribute:
        assert comparison.configuration is config  # matches the object!
        assert comparison in config.comparisons  # back link!

        # Check the query object attribute:
        assert (
            "Example " + hashes[comparison.query_hash] == comparison.query.description
        )

        # Check the subject object attribute:
        assert (
            "Example " + hashes[comparison.subject_hash]
            == comparison.subject.description
        )

    for genome in session.query(db_orm.Genome):
        assert len(genome.query_comparisons) == len(NAMES)
        for comparison in genome.query_comparisons:
            assert comparison.configuration is config
        assert len(genome.subject_comparisons) == len(NAMES)
        for comparison in genome.subject_comparisons:
            assert comparison.configuration is config

    del session  # disconnect

    assert tmp_db.is_file()

    with db_orm.connect_to_db(logger, tmp_db) as new_session:
        assert new_session.query(db_orm.Configuration).count() == 1
        assert new_session.query(db_orm.Genome).count() == len(NAMES)
        assert new_session.query(db_orm.Comparison).count() == len(NAMES) ** 2
        assert new_session.query(db_orm.Run).count() == 0

    tmp_db.unlink()


def test_make_and_populate_runs(tmp_path: str) -> None:
    """Populate new DB with config and runs."""
    tmp_db = Path(tmp_path) / "runs.sqlite"
    assert not tmp_db.is_file()

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)

    config = db_orm.Configuration(
        method="guessing",
        program="guestimate",
        version="v0.1.2beta3",
        fragsize=1000,
        kmersize=31,
    )
    # Test the __repr__
    assert repr(config) == (
        "Configuration(configuration_id=None,"
        " program='guestimate', version='v0.1.2beta3',"
        " fragsize=1000, mode=None, kmersize=31, minmatch=None, extra=None)"
    )
    session.add(config)
    assert config.configuration_id is None
    session.commit()
    assert config.configuration_id == 1

    run_one = db_orm.Run(
        configuration_id=config.configuration_id,
        name="Test One",
        cmdline="pyani_plus run -m guestimate --input ../my-genomes/ -d working.sqlite",
        fasta_directory="../my-genomes/",
        date=datetime.date(2024, 9, 3),
        status="Pending",
    )
    assert repr(run_one) == (
        "Run(run_id=None, configuration_id=1,"
        " cmdline='pyani_plus run -m guestimate --input ../my-genomes/ -d working.sqlite',"
        " date=datetime.date(2024, 9, 3), status='Pending',"
        " name='Test One', ...)"
    )
    session.add(run_one)
    assert run_one.run_id is None
    session.commit()
    assert run_one.run_id == 1

    # These have not been collated yet, and there are in fact no genomes yet either:
    assert run_one.identities is None
    assert run_one.cov_query is None
    assert run_one.aln_length is None
    assert run_one.sim_errors is None
    assert run_one.hadamard is None
    assert run_one.tani is None
    run_one.cache_comparisons()
    # They now exist, but are empty:
    assert run_one.identities is not None
    assert run_one.identities.empty

    run_two = db_orm.Run(
        configuration_id=config.configuration_id,
        name="Test Two",
        cmdline="pyani_plus run -m guestimate --input ../my-genomes/ -d working.sqlite",
        fasta_directory="../my-genomes/",
        date=datetime.date(2024, 9, 4),
        status="Pending",
    )
    assert repr(run_two) == (
        "Run(run_id=None, configuration_id=1,"
        " cmdline='pyani_plus run -m guestimate --input ../my-genomes/ -d working.sqlite',"
        " date=datetime.date(2024, 9, 4), status='Pending',"
        " name='Test Two', ...)"
    )
    session.add(run_two)
    assert run_two.run_id is None
    session.commit()
    assert run_two.run_id == 2  # noqa: PLR2004

    assert run_one.configuration is run_two.configuration
    # Now check can access all the runs from a configuration object
    runs = list(config.runs)
    assert len(runs) == 2  # noqa: PLR2004
    assert runs[0] is run_one
    assert runs[1] is run_two

    del session
    assert tmp_db.is_file()
    with db_orm.connect_to_db(logger, tmp_db) as new_session:
        assert new_session.query(db_orm.Configuration).count() == 1
        assert new_session.query(db_orm.Genome).count() == 0
        assert new_session.query(db_orm.Comparison).count() == 0
        assert new_session.query(db_orm.Run).count() == 2  # noqa: PLR2004
    tmp_db.unlink()


def test_make_and_populate_mock_example(tmp_path: str) -> None:  # noqa: PLR0915
    """Populate new DB with config, runs, genomes and comparisons."""
    tmp_db = Path(tmp_path) / "mock.sqlite"
    assert not tmp_db.is_file()

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)

    config = db_orm.Configuration(
        method="guessing",
        program="guestimate",
        version="v0.1.2beta3",
        fragsize=1000,
        kmersize=31,
    )
    assert repr(config) == (
        "Configuration(configuration_id=None,"
        " program='guestimate', version='v0.1.2beta3',"
        " fragsize=1000, mode=None, kmersize=31, minmatch=None, extra=None)"
    )
    session.add(config)
    session.commit()

    run = db_orm.Run(
        configuration_id=config.configuration_id,
        name="Empty",
        cmdline="pyani_plus run -m guestimate --input ../my-genomes/ -d working.sqlite",
        fasta_directory="../my-genomes/",
        date=datetime.date(2023, 12, 25),
        status="Aborted",
    )
    assert repr(run) == (
        "Run(run_id=None, configuration_id=1,"
        " cmdline='pyani_plus run -m guestimate --input ../my-genomes/ -d working.sqlite',"
        " date=datetime.date(2023, 12, 25), status='Aborted',"
        " name='Empty', ...)"
    )
    session.add(run)

    run = db_orm.Run(
        configuration_id=config.configuration_id,
        name="Test Run",
        cmdline="pyani_plus run -m guestimate --input ../my-genomes/ -d working.sqlite",
        fasta_directory="../my-genomes/",
        date=datetime.date(2023, 12, 25),
        status="Complete",
    )
    assert repr(run) == (
        "Run(run_id=None, configuration_id=1,"
        " cmdline='pyani_plus run -m guestimate --input ../my-genomes/ -d working.sqlite',"
        " date=datetime.date(2023, 12, 25), status='Complete',"
        " name='Test Run', ...)"
    )
    session.add(run)

    hashes = []
    # Going to record 4 genomes, and all 4x4=16 comparisons
    # However, only going to link 2 genomes to the run (so 4 comparisons)
    for name in ("Genome A", "Genome C", "Genome G", "Genome T"):
        seq = (name[-1] + "ACGT") * 1000
        md5 = str_md5sum(f">{name}\n{seq}\n")
        hashes.append(md5)
        genome = db_orm.Genome(
            genome_hash=md5,
            path=f"../my-genomes/{name}.fasta",
            length=len(seq),
            description=name,
        )
        assert repr(genome) == (
            f"Genome(genome_hash={md5!r}, path='../my-genomes/{name}.fasta',"
            f" length={len(seq)}, description='{name}')"
        )
        session.add(genome)
        if name[-1] in ("A", "T"):
            # Don't know the run_id until the run is committed
            session.add(
                db_orm.RunGenomeAssociation(
                    run=run,
                    genome_hash=md5,
                    fasta_filename=f"{name}.fasta",  # no path here, just filename
                )
            )

    for a in hashes:
        for b in hashes:
            comparison = db_orm.Comparison(
                configuration_id=config.configuration_id,
                query_hash=a,
                subject_hash=b,
                identity=0.99 if a == b else 0.96,
                aln_length=4996 if a == b else 4975,
                uname_system="Darwin",
                uname_release="21.6.0",
                uname_machine="arm64",
            )
            assert repr(comparison) == (
                "Comparison(comparison_id=None,"
                f" query_hash={a!r}, subject_hash={b!r},"
                f" configuration_id={config.configuration_id},"
                f" identity={0.99 if a == b else 0.96},"
                f" aln_length={4996 if a == b else 4975},"
                " sim_errors=None, cov_query=None, cov_subject=None,"
                " uname_system='Darwin', uname_release='21.6.0', uname_machine='arm64')"
            )
            session.add(comparison)

    # These have not been collated yet:
    assert run.identities is None
    run.cache_comparisons()
    assert run.identities is not None

    # The run has only 2 genomes in it by construction
    assert run.identities.equals(
        pd.DataFrame(
            data=np.array([[0.99, 0.96], [0.96, 0.99]], float),
            index=[
                "0ac5f24c37ea8ef2d0e37fbfa61e2a43",
                "4bce437e7bdd91e35d18bfe294dee207",
            ],
            columns=[
                "0ac5f24c37ea8ef2d0e37fbfa61e2a43",
                "4bce437e7bdd91e35d18bfe294dee207",
            ],
        )
    )
    assert run.cov_query.equals(
        pd.DataFrame(
            data=np.array([[np.nan, np.nan], [np.nan, np.nan]], float),
            index=[
                "0ac5f24c37ea8ef2d0e37fbfa61e2a43",
                "4bce437e7bdd91e35d18bfe294dee207",
            ],
            columns=[
                "0ac5f24c37ea8ef2d0e37fbfa61e2a43",
                "4bce437e7bdd91e35d18bfe294dee207",
            ],
        )
    )
    assert run.aln_length.equals(
        pd.DataFrame(
            data=np.array([[4996, 4975], [4975, 4996]], int),
            index=[
                "0ac5f24c37ea8ef2d0e37fbfa61e2a43",
                "4bce437e7bdd91e35d18bfe294dee207",
            ],
            columns=[
                "0ac5f24c37ea8ef2d0e37fbfa61e2a43",
                "4bce437e7bdd91e35d18bfe294dee207",
            ],
        )
    )
    assert run.sim_errors.equals(
        pd.DataFrame(
            data=np.array([[np.nan, np.nan], [np.nan, np.nan]], float),
            index=[
                "0ac5f24c37ea8ef2d0e37fbfa61e2a43",
                "4bce437e7bdd91e35d18bfe294dee207",
            ],
            columns=[
                "0ac5f24c37ea8ef2d0e37fbfa61e2a43",
                "4bce437e7bdd91e35d18bfe294dee207",
            ],
        )
    )
    # float * nan = nan, so hadamard is all nan:
    assert run.hadamard.equals(
        pd.DataFrame(
            data=np.array([[np.nan, np.nan], [np.nan, np.nan]], float),
            index=[
                "0ac5f24c37ea8ef2d0e37fbfa61e2a43",
                "4bce437e7bdd91e35d18bfe294dee207",
            ],
            columns=[
                "0ac5f24c37ea8ef2d0e37fbfa61e2a43",
                "4bce437e7bdd91e35d18bfe294dee207",
            ],
        )
    )

    session.commit()

    # For debug testing with sqlite3: import shutil; shutil.copy(tmp_db, "demo.sqlite")
    del session, config, run, genome, comparison
    assert tmp_db.is_file()
    with db_orm.connect_to_db(logger, tmp_db) as new_session:
        config = new_session.query(db_orm.Configuration).one()
        assert new_session.query(db_orm.Genome).count() == 4  # noqa: PLR2004
        genomes = set(new_session.query(db_orm.Genome))
        assert new_session.query(db_orm.Comparison).count() == 4 * 4
        run = new_session.query(db_orm.Run).where(db_orm.Run.status == "Complete").one()
        assert run.configuration is config
        assert set(run.genomes).issubset(genomes)  # order seemed to be stochastic
        for comparison in new_session.query(db_orm.Comparison):
            assert comparison.configuration is config
            assert comparison.query in genomes
            assert comparison.subject in genomes
        # Confirm we can pull out only the comparisons belonging to a run:
        assert len(list(run.comparisons())) == 2 * 2  # only 4, not all 16
        #
        # Using session with echo=True (and making the test fail to see the output),
        # could confirm the SQL matches my goal:
        #
        # SELECT ... FROM comparisons
        # JOIN runs_genomes AS run_query ON comparisons.query_hash = run_query.genome_hash
        # JOIN runs_genomes AS run_subject ON comparisons.subject_hash = run_subject.genome_hash
        # WHERE run_query.run_id = ? AND run_subject.run_id = ?
    tmp_db.unlink()


def test_add_config(tmp_path: str) -> None:
    """Confirm repeating adding a configuration has no effect."""
    tmp_db = Path(tmp_path) / "config.sqlite"
    assert not tmp_db.is_file()

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)

    with pytest.raises(
        NoResultFound,
        match="Requested configuration not already in DB",
    ):
        db_orm.db_configuration(
            session,
            method="guessing",
            program="guestimate",
            version="v0.1.2beta3",
            fragsize=100,
            kmersize=17,
            create=False,
        )

    config = db_orm.db_configuration(
        session,
        method="guessing",
        program="guestimate",
        version="v0.1.2beta3",
        fragsize=100,
        kmersize=17,
        create=True,
    )
    assert repr(config) == (
        "Configuration(configuration_id=1,"
        " program='guestimate', version='v0.1.2beta3',"
        " fragsize=100, mode=None, kmersize=17, minmatch=None, extra=None)"
    )

    # Trying to add the exact same values should return the existing entry:
    assert config is db_orm.db_configuration(
        session,
        method="guessing",
        program="guestimate",
        version="v0.1.2beta3",
        fragsize=100,
        kmersize=17,
        create=True,
    )


def test_add_genome(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Confirm repeating adding a genome has no effect."""
    tmp_db = Path(tmp_path) / "config.sqlite"
    assert not tmp_db.is_file()

    fasta = next(input_genomes_tiny.glob("*.f*"))

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    md5 = file_md5sum(fasta)
    with pytest.raises(
        NoResultFound,
        match="Requested genome not already in DB",
    ):
        db_orm.db_genome(logger, session, fasta, md5, create=False)

    genome = db_orm.db_genome(logger, session, fasta, md5, create=True)
    # Adding again should return the original row again
    assert genome is db_orm.db_genome(logger, session, fasta, md5, create=True)


def test_add_genome_not_gzipped(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Confirm catches an uncompressed FASTA with .gz extension."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "not_gzip.sqlite"
    tmp_input = tmp_dir / "input"
    tmp_input.mkdir()

    file = next(input_genomes_tiny.glob("*.f*"))
    fasta = tmp_input / (file.name + ".gz")
    fasta.symlink_to(file)

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    md5 = file_md5sum(fasta)
    with pytest.raises(
        SystemExit,
        match="Has .gz ending, but .*\\.gz is NOT gzip compressed",
    ):
        db_orm.db_genome(logger, session, fasta, md5, create=True)


def test_add_genome_no_gz_ext(tmp_path: str, input_gzip_bacteria: Path) -> None:
    """Confirm catches a compressed FASTA without .gz extension."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "no_gz.sqlite"
    tmp_input = tmp_dir / "input"
    tmp_input.mkdir()

    file = next(input_gzip_bacteria.glob("*.f*.gz"))
    fasta = tmp_input / (file.name[:-3])
    fasta.symlink_to(file)

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    md5 = file_md5sum(fasta)
    with pytest.raises(
        SystemExit,
        match="No .gz ending, but .*\\.f.* is gzip compressed",
    ):
        db_orm.db_genome(logger, session, fasta, md5, create=True)


def test_helper_functions(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Populate new DB using helper functions."""
    tmp_db = Path(tmp_path) / "mock.sqlite"
    assert not tmp_db.is_file()

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)

    config = db_orm.db_configuration(
        session,
        method="guessing",
        program="guestimate",
        version="v0.1.2beta3",
        fragsize=1000,
        kmersize=31,
        create=True,
    )
    assert repr(config) == (
        "Configuration(configuration_id=1,"
        " program='guestimate', version='v0.1.2beta3',"
        " fragsize=1000, mode=None, kmersize=31, minmatch=None, extra=None)"
    )
    fasta_to_hash = {}
    for fasta_filename in input_genomes_tiny.glob("*.f*"):
        md5 = file_md5sum(fasta_filename)
        db_orm.db_genome(logger, session, fasta_filename, md5, create=True)
        fasta_to_hash[fasta_filename] = md5

    run = db_orm.add_run(
        session,
        configuration=config,
        fasta_to_hash=fasta_to_hash,
        status="Started",
        fasta_directory=Path("/mnt/shared/genomes/"),
        name="Guess Run",
        date=None,
        cmdline="pyani_plus run --method guestimate --fasta blah blah",
    )
    now = run.date
    assert repr(run) == (
        "Run(run_id=1, configuration_id=1, "
        "cmdline='pyani_plus run --method guestimate --fasta blah blah', "
        f"date={now!r}, status='Started', name='Guess Run', ...)"
    )
    assert run.genomes.count() == len(fasta_to_hash)
    assert run.comparisons().count() == 0  # not logged yet

    # At this point in a real run we would start parallel worker jobs
    # to compute the comparisons (possible spread over a cluster):
    for a in fasta_to_hash.values():
        for b in fasta_to_hash.values():
            db_orm.db_comparison(
                session,
                config.configuration_id,
                a,
                b,
                1 if a == b else 0.99,
                12345,
                cov_query=0.99,
                cov_subject=0.99,
            )

    # Now that all the comparisons are in the DB, can collate and cache matrices
    assert run.comparisons().count() == len(fasta_to_hash) ** 2
    run.cache_comparisons()
    run.status = "Complete"
    session.commit()
    assert run.df_hadamard == (
        '{"columns":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],"index":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],"data":[[0.99,0.9801,0.9801],[0.9801,0.99,0.9801],[0.9801,0.9801,0.99]]}'
    )
    tmp_db.unlink()


def test_insert_no_comps(tmp_path: str) -> None:
    """Checking a corner case with recording no comparisons."""
    tmp_db = Path(tmp_path) / "empty.db"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    assert db_orm.insert_comparisons_with_retries(logger, session, [], "test only")
    session.close()
    # Is there an easy way to test if this called commit or not?
