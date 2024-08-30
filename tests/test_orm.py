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
    names = ("Alpha", "Beta", "Gamma")
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

    for name in names:
        seq = "ACGT" * 100 * len(name)
        fasta = f">{name}\n{seq}\n"
        md5 = hashlib.md5(fasta.encode("ascii")).hexdigest()  # noqa: S324
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
    assert session.query(db_orm.Genome).count() == len(names)
    del session  # disconnect

    assert tmp_db.is_file()
    with tmp_db.open("rb") as handle:
        magic = handle.read(16)
        assert magic == b"SQLite format 3\0"

    with db_orm.connect_to_db(tmp_db) as new_session:
        assert new_session.query(db_orm.Configuration).count() == 1
        assert new_session.query(db_orm.Genome).count() == len(names)

    tmp_db.unlink()
