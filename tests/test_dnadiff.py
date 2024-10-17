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
"""Test methods for calculating dnadiff.

These tests are intended to be run from the repository root using:

make test
"""

# Required to support pytest automated testing
from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, tools
from pyani_plus.methods import method_dnadiff

from . import get_matrix_entry


@pytest.fixture
def expected_mcoords_output() -> tuple[float, int]:
    """Expected average identity and aligned bases with gaps.

    MGV-GENOME-0264574 (REF) MGV-GENOME-0266457 (QRY).
    """  # noqa: D401
    return (99.63, 39594)


@pytest.fixture
def expected_gap_lengths_qry() -> int:
    """Expected gap length in the QRY alignment."""  # noqa: D401
    return 418


def test_parse_mcoords(
    dnadiff_targets_showcoords_indir: Path, expected_mcoords_output: tuple[float, int]
) -> None:
    """Check parsing of test mcoords file."""
    assert expected_mcoords_output == method_dnadiff.parse_mcoords(
        dnadiff_targets_showcoords_indir
        / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.mcoords"
    )


def test_parse_qdiff(
    dnadiff_targets_showdiff_indir: Path, expected_gap_lengths_qry: int
) -> None:
    """Check parsing of test qdiff file."""
    assert expected_gap_lengths_qry == method_dnadiff.parse_qdiff(
        dnadiff_targets_showdiff_indir
        / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.qdiff"
    )


def test_missing_db(
    tmp_path: str,
    input_genomes_tiny: Path,
    dnadiff_targets_showcoords_indir: Path,
    dnadiff_targets_showdiff_indir: Path,
) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.log_dnadiff(
            database=tmp_db,
            # These are for the comparison table
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            mcoords=dnadiff_targets_showcoords_indir
            / "MGV-GENOME-0266457_vs_MGV-GENOME-0264574.mcoords",
            qdiff=dnadiff_targets_showdiff_indir
            / "MGV-GENOME-0266457_vs_MGV-GENOME-0264574.qdiff",
        )


def test_bad_query_or_subject(
    tmp_path: str,
    input_genomes_tiny: Path,
    dnadiff_targets_showcoords_indir: Path,
    dnadiff_targets_showdiff_indir: Path,
) -> None:
    """Mismatch between query or subject FASTA in show-coords and show-diff output and commandline."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --query-fasta .*/MGV-GENOME-0266457.fna"
            " but query in mcoords filename was MGV-GENOME-0264574"
        ),
    ):
        private_cli.log_dnadiff(
            database=tmp_db,
            query_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            mcoords=dnadiff_targets_showcoords_indir
            / "MGV-GENOME-0266457_vs_MGV-GENOME-0264574.mcoords",
            qdiff=dnadiff_targets_showdiff_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.qdiff",
        )

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --subject-fasta .*/MGV-GENOME-0264574.fas"
            " but query in mcoords filename was MGV-GENOME-0266457"
        ),
    ):
        private_cli.log_dnadiff(
            database=tmp_db,
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            mcoords=dnadiff_targets_showcoords_indir
            / "MGV-GENOME-0266457_vs_MGV-GENOME-0264574.mcoords",
            qdiff=dnadiff_targets_showdiff_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.qdiff",
        )

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --query-fasta .*/OP073605.fasta"
            " but query in qdiff filename was MGV-GENOME-0266457"
        ),
    ):
        private_cli.log_dnadiff(
            database=tmp_db,
            query_fasta=input_genomes_tiny / "OP073605.fasta",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            mcoords=dnadiff_targets_showcoords_indir
            / "MGV-GENOME-0266457_vs_OP073605.mcoords",
            qdiff=dnadiff_targets_showdiff_indir
            / "OP073605_vs_MGV-GENOME-0266457.qdiff",
        )

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --subject-fasta .*/MGV-GENOME-0266457.fna"
            " but query in qdiff filename was MGV-GENOME-0264574"
        ),
    ):
        private_cli.log_dnadiff(
            database=tmp_db,
            query_fasta=input_genomes_tiny / "OP073605.fasta",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            mcoords=dnadiff_targets_showcoords_indir
            / "MGV-GENOME-0266457_vs_OP073605.mcoords",
            qdiff=dnadiff_targets_showdiff_indir
            / "MGV-GENOME-0264574_vs_OP073605.qdiff",
        )


def test_logging_dnadiff(
    tmp_path: str,
    input_genomes_tiny: Path,
    dnadiff_targets_showcoords_indir: Path,
    dnadiff_targets_showdiff_indir: Path,
    dir_dnadiff_matrices: Path,
) -> None:
    """Check can log a fastANI comparison to DB."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_dnadiff()

    private_cli.log_configuration(
        database=tmp_db,
        method="dnadiff",
        program=tool.exe_path.stem,
        version=tool.version,
        create_db=True,
    )
    private_cli.log_genome(
        database=tmp_db,
        fasta=[
            input_genomes_tiny / "MGV-GENOME-0264574.fas",
            input_genomes_tiny / "MGV-GENOME-0266457.fna",
        ],
    )

    private_cli.log_dnadiff(
        database=tmp_db,
        query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
        subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
        mcoords=dnadiff_targets_showcoords_indir
        / "MGV-GENOME-0266457_vs_MGV-GENOME-0264574.mcoords",
        qdiff=dnadiff_targets_showdiff_indir
        / "MGV-GENOME-0266457_vs_MGV-GENOME-0264574.qdiff",
    )

    # Check the recorded comparison values
    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 1
    comp = session.query(db_orm.Comparison).one()
    query = "689d3fd6881db36b5e08329cf23cecdd"  # MGV-GENOME-0264574.fas
    subject = "78975d5144a1cd12e98898d573cf6536"  # MGV-GENOME-0266457.fna
    pytest.approx(
        comp.identity,
        get_matrix_entry(
            dir_dnadiff_matrices / "avg_identity_matrix.tsv", query, subject
        ),
    )
    pytest.approx(
        comp.aln_length,
        get_matrix_entry(
            dir_dnadiff_matrices / "aligned_bases_matrix.tsv", query, subject
        ),
    )
    pytest.approx(
        comp.cov_query,
        get_matrix_entry(dir_dnadiff_matrices / "coverage_matrix.tsv", query, subject),
    )
    session.close()
    tmp_db.unlink()
