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
"""Tests for the sourmash implementation.

These tests are intended to be run from the repository root using:

make test
"""

# Required to support pytest automated testing
from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, tools
from pyani_plus.methods import method_sourmash

from . import get_matrix_entry


def test_compare_parsing(sourmash_targets_compare_indir: Path) -> None:
    """Check parsing of sourmash compare .csv file."""
    assert (
        method_sourmash.parse_compare(
            sourmash_targets_compare_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv"
        )
        == 0.997900938305757  # noqa: PLR2004
    )


def test_missing_db(
    tmp_path: str, input_genomes_tiny: Path, sourmash_targets_compare_indir: Path
) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.log_sourmash(
            database=tmp_db,
            # These are for the comparison table
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            compare=sourmash_targets_compare_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv",
        )


def test_bad_query_or_subject(
    tmp_path: str, input_genomes_tiny: Path, sourmash_targets_compare_indir: Path
) -> None:
    """Mismatch between query or subject FASTA in sourmash compare output and commandline."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --query-fasta .*/MGV-GENOME-0266457.fna"
            " but query in sourmash compare filename was MGV-GENOME-0264574"
        ),
    ):
        private_cli.log_sourmash(
            database=tmp_db,
            query_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            compare=sourmash_targets_compare_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv",
        )

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --subject-fasta .*/MGV-GENOME-0264574.fas"
            " but subject in sourmash compare filename was MGV-GENOME-0266457"
        ),
    ):
        private_cli.log_sourmash(
            database=tmp_db,
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            compare=sourmash_targets_compare_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv",
        )


def test_logging_sourmash(
    tmp_path: str,
    input_genomes_tiny: Path,
    sourmash_targets_compare_indir: Path,
    dir_sourmash_results: Path,
) -> None:
    """Check can log a sourmash comparison to DB."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_sourmash()

    private_cli.log_configuration(
        database=tmp_db,
        method="sourmash",
        program=tool.exe_path.stem,
        version=tool.version,
        mode=method_sourmash.MODE,
        kmersize=method_sourmash.KMER_SIZE,
        extra="scaled=" + str(method_sourmash.SCALED),
        create_db=True,
    )
    private_cli.log_genome(
        database=tmp_db,
        fasta=[
            input_genomes_tiny / "MGV-GENOME-0264574.fas",
            input_genomes_tiny / "MGV-GENOME-0266457.fna",
        ],
    )

    private_cli.log_sourmash(
        database=tmp_db,
        query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
        subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
        compare=sourmash_targets_compare_indir
        / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv",
    )

    # Check the recorded comparison values
    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 1
    comp = session.query(db_orm.Comparison).one()
    pytest.approx(
        comp.identity,
        get_matrix_entry(
            dir_sourmash_results / "matrix_identity.tsv",
            comp.query_hash,
            comp.subject_hash,
        ),
    )
    session.close()
    tmp_db.unlink()
