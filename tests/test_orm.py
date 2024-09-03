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
"""Tests for the pyani_plus/db_orm.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

import datetime
import hashlib
import platform
from pathlib import Path

from pyani_plus import db_orm


def test_make_new_db(tmp_path: str) -> None:
    """Confirm can create a new empty database."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()
    db_orm.connect_to_db(tmp_db)  # discard the session, should close
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

    session = db_orm.connect_to_db(tmp_db)

    uname = platform.uname()
    config = db_orm.Configuration(
        method="guessing",
        machine=uname.machine,  # CPU arch
        system=uname.system,  # Operating system
        program="guestimate",
        version="v0.1.2beta3",
        fragsize=100,
        kmersize=17,
    )
    # Test the __repr__
    assert repr(config) == (
        "Configuration(configuration_id=None,"
        f" machine={uname.machine!r}, system={uname.system!r},"
        " program='guestimate', version='v0.1.2beta3',"
        " fragsize=100, maxmatch=None, kmersize=17, minmatch=None)"
    )
    session.add(config)
    assert config.configuration_id is None
    session.commit()
    assert config.configuration_id == 1

    hashes = {}
    for name in NAMES:
        seq = "ACGT" * int(DUMMY_ALIGN_LEN / 4) * len(name)
        fasta = f">{name}\n{seq}\n"
        md5 = hashlib.md5(fasta.encode("ascii")).hexdigest()  # noqa: S324
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
            )
            assert comparison.configuration_id == config.configuration_id
            assert repr(comparison) == (
                "Comparison(comparison_id=None,"
                f" query_hash={query!r}, subject_hash={subject!r},"
                f" configuration_id={config.configuration_id},"
                f" identity=0.96, aln_length={DUMMY_ALIGN_LEN}, sim_errs=None,"
                " cov_query=None, cov_subject=None)"
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
            "FragSize: 100, MaxMatch: None, KmerSize: 17, MinMatch: None"
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

    with db_orm.connect_to_db(tmp_db) as new_session:
        assert new_session.query(db_orm.Configuration).count() == 1
        assert new_session.query(db_orm.Genome).count() == len(NAMES)
        assert new_session.query(db_orm.Comparison).count() == len(NAMES) ** 2
        assert new_session.query(db_orm.Run).count() == 0

    tmp_db.unlink()


def test_make_and_populate_runs(tmp_path: str) -> None:
    """Populate new DB with config and runs."""
    tmp_db = Path(tmp_path) / "runs.sqlite"
    assert not tmp_db.is_file()

    session = db_orm.connect_to_db(tmp_db)

    uname = platform.uname()
    config = db_orm.Configuration(
        method="guessing",
        machine=uname.machine,  # CPU arch
        system=uname.system,  # Operating system
        program="guestimate",
        version="v0.1.2beta3",
        fragsize=1000,
        kmersize=31,
    )
    # Test the __repr__
    assert repr(config) == (
        "Configuration(configuration_id=None,"
        f" machine={uname.machine!r}, system={uname.system!r},"
        " program='guestimate', version='v0.1.2beta3',"
        " fragsize=1000, maxmatch=None, kmersize=31, minmatch=None)"
    )
    session.add(config)
    assert config.configuration_id is None
    session.commit()
    assert config.configuration_id == 1

    run_one = db_orm.Run(
        configuration_id=config.configuration_id,
        name="Test One",
        cmdline="pyani_plus run -m guestimate --input ../my-genomes/ -d working.sqlite",
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

    run_two = db_orm.Run(
        configuration_id=config.configuration_id,
        name="Test Two",
        cmdline="pyani_plus run -m guestimate --input ../my-genomes/ -d working.sqlite",
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
    with db_orm.connect_to_db(tmp_db) as new_session:
        assert new_session.query(db_orm.Configuration).count() == 1
        assert new_session.query(db_orm.Genome).count() == 0
        assert new_session.query(db_orm.Comparison).count() == 0
        assert new_session.query(db_orm.Run).count() == 2  # noqa: PLR2004
    tmp_db.unlink()
