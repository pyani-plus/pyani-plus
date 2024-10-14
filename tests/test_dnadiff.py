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

import pandas as pd
import pytest

from pyani_plus import utils
from pyani_plus.methods import method_dnadiff


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


@pytest.fixture
def expected_dnadiff_results() -> method_dnadiff.ComparisonResultDnadiff:
    """Example of expected aniM ComparisionResult."""  # noqa: D401
    return method_dnadiff.ComparisonResultDnadiff(
        qname="MGV-GENOME-0266457",
        rname="MGV-GENOME-0264574",
        q_aligned_bases_with_gaps=39594,
        avg_id=99.63,
        q_length=39594,
        r_length=39253,
        alignment_gaps=418,
        aligned_bases=39176,
        q_cov=(39176 / 39594) * 100,
    )


def compare_results(
    dataframe: pd.DataFrame,
    mcoords: Path,
    dnadiff_targets_showdiff_indir: Path,
    input_genomes_tiny: Path,
    item: str,
) -> bool:
    """Compare dnadiff results.

    Compare values as they are, except for coverage values,
    which are rounded to 2 decimal places.
    """
    genomes = list(input_genomes_tiny.glob("*.f*"))
    genome_hashes = {record.stem: utils.file_md5sum(record) for record in genomes}
    reference, query = mcoords.stem.split("_vs_")

    result = method_dnadiff.collect_dnadiff_results(
        mcoords,
        dnadiff_targets_showdiff_indir / f"{reference}_vs_{query}.qdiff",
        input_genomes_tiny,
    )

    if item == "q_cov":
        result_value = round(result.item(item), 2)
        dataframe_value = round(
            dataframe.loc[genome_hashes[reference], genome_hashes[query]], 2
        )
    else:
        result_value = result.item(item)
        dataframe_value = dataframe.loc[genome_hashes[reference], genome_hashes[query]]

    return result_value == dataframe_value


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


def test_collect_results(
    expected_dnadiff_results: method_dnadiff.ComparisonResultDnadiff,
    dnadiff_targets_showdiff_indir: Path,
    dnadiff_targets_showcoords_indir: Path,
    input_genomes_tiny: Path,
) -> None:
    """Check that ComparisionResultDnadiff class is populated with correct values."""
    assert expected_dnadiff_results == method_dnadiff.collect_dnadiff_results(
        dnadiff_targets_showcoords_indir
        / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.mcoords",
        dnadiff_targets_showdiff_indir
        / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.qdiff",
        input_genomes_tiny,
    )


def test_aligned_bases(
    dir_dnadiff_matrices: Path,
    dnadiff_targets_showcoords_indir: Path,
    input_genomes_tiny: Path,
    dnadiff_targets_showdiff_indir: Path,
) -> None:
    """Test dnadiff AlignedBases values."""
    for fname in (dnadiff_targets_showcoords_indir).glob("*.mcoords"):
        assert compare_results(
            pd.read_csv(
                dir_dnadiff_matrices / "aligned_bases_matrix.tsv", sep="\t", index_col=0
            ),
            fname,
            dnadiff_targets_showdiff_indir,
            input_genomes_tiny,
            "aligned_bases",
        )


def test_avg_identity(
    dir_dnadiff_matrices: Path,
    dnadiff_targets_showcoords_indir: Path,
    input_genomes_tiny: Path,
    dnadiff_targets_showdiff_indir: Path,
) -> None:
    """Test dnadiff average identity values."""
    for fname in (dnadiff_targets_showcoords_indir).glob("*.mcoords"):
        assert compare_results(
            pd.read_csv(
                dir_dnadiff_matrices / "avg_identity_matrix.tsv", sep="\t", index_col=0
            ),
            fname,
            dnadiff_targets_showdiff_indir,
            input_genomes_tiny,
            "avg_id",
        )


def test_coverage(
    dir_dnadiff_matrices: Path,
    dnadiff_targets_showcoords_indir: Path,
    input_genomes_tiny: Path,
    dnadiff_targets_showdiff_indir: Path,
) -> None:
    """Test dnadiff average identity values."""
    for fname in (dnadiff_targets_showcoords_indir).glob("*.mcoords"):
        assert compare_results(
            pd.read_csv(
                dir_dnadiff_matrices / "coverage_matrix.tsv", sep="\t", index_col=0
            ),
            fname,
            dnadiff_targets_showdiff_indir,
            input_genomes_tiny,
            "q_cov",
        )
