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
"""Tests for the pyani_plus/db_orm.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

import hashlib
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


def test_make_and_populate_new_db(tmp_path: str) -> None:
    """Confirm can create a new empty database, add genomes, close, and reopen."""
    DUMMY_ALIGN_LEN = 400  # noqa: N806
    NAMES = ("Alpha", "Beta", "Gamma")  # noqa: N806
    tmp_db = Path(tmp_path) / "genomes.sqlite"
    assert not tmp_db.is_file()

    session = db_orm.connect_to_db(tmp_db)

    config = db_orm.Configuration(
        method="guessing",
        program="guestimate",
        version="v0.1.2beta3",
        fragsize=17,
        kmersize=17,
    )
    # Test the __repr__
    assert repr(config) == (
        "Configuration(configuration_id=None,"
        " program='guestimate', version='v0.1.2beta3',"
        " fragsize=17, maxmatch=None, kmersize=17, minmatch=None)"
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
    session.commit()
    assert len(config.comparisons) == 0
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
                f"Comparison(comparison_id=None,"
                f" query_hash={query!r}, subject_hash={subject!r},"
                f" configuration_id={config.configuration_id},"
                f" identity=0.96, aln_length={DUMMY_ALIGN_LEN}, sim_errs=None,"
                " cov_query=None, cov_subject=None)"
            )
            session.add(comparison)
    session.commit()
    assert session.query(db_orm.Comparison).count() == len(NAMES) ** 2
    assert len(config.comparisons) == len(NAMES) ** 2
    for comparison in config.comparisons:
        assert comparison.aln_length == DUMMY_ALIGN_LEN

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
    with tmp_db.open("rb") as handle:
        magic = handle.read(16)
        assert magic == b"SQLite format 3\0"

    with db_orm.connect_to_db(tmp_db) as new_session:
        assert new_session.query(db_orm.Configuration).count() == 1
        assert new_session.query(db_orm.Genome).count() == len(NAMES)
        assert new_session.query(db_orm.Comparison).count() == len(NAMES) ** 2

    tmp_db.unlink()
