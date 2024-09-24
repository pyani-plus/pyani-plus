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

from pathlib import Path

import pytest

from pyani_plus import private_cli


def test_log_configuration(tmp_path: str) -> None:
    """Confirm can create a new empty database via log-configuration."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit) as err:
        private_cli.log_configuration(
            tmp_db,
            method="guessing",
            program="guestimate",
            version="0.1.2beta3",
            fragsize=100,
            kmersize=51,
            create_db=False,
        )
    assert "does not exist, but not using --create-db" in str(err)

    # This time create it
    private_cli.log_configuration(
        tmp_db,
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        create_db=True,
    )

    # This time should already be a DB there
    private_cli.log_configuration(
        tmp_db,
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=75,
        kmersize=31,
        create_db=False,
    )

    tmp_db.unlink()


def test_log_genome(tmp_path: str, input_genomes_small: Path) -> None:
    """Confirm can create a new empty database via log-genome."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit) as err:
        private_cli.log_genome(
            database=tmp_db,
            fasta=list(
                input_genomes_small.glob("*.fasta")  # subset of folder
            ),
        )
    assert "does not exist, but not using --create-db" in str(err)

    # This time create it
    private_cli.log_genome(
        database=tmp_db,
        fasta=list(
            input_genomes_small.glob("*.fasta")  # subset of folder
        ),
        create_db=True,
    )


def test_log_comparison(tmp_path: str, input_genomes_small: Path) -> None:
    """Confirm can create a mock DB using log-comparison etc."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    private_cli.log_configuration(
        tmp_db,
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        create_db=True,
    )

    fasta = list(input_genomes_small.glob("*.fasta"))  # subset of folder
    private_cli.log_genome(
        database=tmp_db,
        fasta=fasta,
        create_db=False,
    )

    # Should at this point log the run

    for query in fasta:
        for subject in fasta:
            private_cli.log_comparison(
                database=tmp_db,
                query_fasta=query,
                subject_fasta=subject,
                identity=1.0 if query == subject else 0.96,
                aln_length=12345,
                method="guessing",
                program="guestimate",
                version="0.1.2beta3",
                fragsize=100,
                kmersize=51,
                sim_errors=1,
                cov_query=0.98,
                cov_subject=0.98,
                create_db=False,
            )

    # Would next test updating the run matrices
