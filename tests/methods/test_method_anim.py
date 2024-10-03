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
"""Test methods for calculating ANIm.

These tests are intended to be run from the repository root using:

make test
"""

# Required to support pytest automated testing
from pathlib import Path

import pandas as pd
import pytest

from pyani_plus.methods import method_anim


@pytest.fixture
def aligned_regions() -> dict:
    """Example of aligned regions with overlaps."""  # noqa: D401
    return {
        "MGV-GENOME-0266457": [(1, 37636), (17626, 39176)],
    }


@pytest.fixture
def expected_anim_results() -> method_anim.ComparisonResult:
    """Example of expected aniM ComparisionResult."""  # noqa: D401
    return method_anim.ComparisonResult(
        qname="MGV-GENOME-0264574",
        rname="MGV-GENOME-0266457",
        q_aligned_bases=39169,
        r_aligned_bases=39176,
        sim_errs=222,
        avg_id=0.9962487643734,
        program="nucmer",
        q_length=39253,
        r_length=39594,
        q_cov=39169 / 39253,
        r_cov=39176 / 39594,
        r_hadamard=(39176 / 39594) * 0.9962487643734,
        q_hadamard=(39169 / 39253) * 0.9962487643734,
    )


def compare_results_decimal(filterfile: Path, datadir: pd.DataFrame, val: int) -> bool:
    """Compare delta-filter files and expected anim results.

    This function expects delta-filter file and dataframe with expected
    anim results as input and returns True if the anim values is the same,
    and False if the two files differ. It also expects index specifying
    which value from the parse_delta function to compare, as the parse_delta
    function returns several values, including aligned bases in reference,
    aligned bases in query, aligned identity, and similarity errors.
    """
    reference, query = filterfile.stem.split("_vs_")
    pyani_label = {_.split(":")[0]: _ for _ in datadir}

    return round(method_anim.parse_delta(filterfile)[val], 6) == round(
        datadir.loc[pyani_label[reference], pyani_label[query]], 6
    )


def compare_results_no_decimal(
    filterfile: Path, datadir: pd.DataFrame, val: int
) -> bool:
    """Compare delta-filter files and expected anim results.

    Without rounding numbers to 6 decimal places.
    """
    reference, query = filterfile.stem.split("_vs_")
    pyani_label = {_.split(":")[0]: _ for _ in datadir}

    return (
        method_anim.parse_delta(filterfile)[val]
        == datadir.loc[pyani_label[reference], pyani_label[query]]
    )


def test_average_identity(
    anim_nucmer_targets_filter_indir: Path,
    dir_anim_results: Path,
) -> None:
    """Check aniM average identity."""
    for fname in (anim_nucmer_targets_filter_indir).glob("*.filter"):
        assert compare_results_decimal(
            fname,
            pd.read_csv(
                dir_anim_results / "matrix_identity_1.tab", index_col=0, sep="\t"
            ),
            2,
        )


def test_aln_lengths(
    anim_nucmer_targets_filter_indir: Path,
    dir_anim_results: Path,
) -> None:
    """Check aniM alignment lengths."""
    for fname in (anim_nucmer_targets_filter_indir).glob("*.filter"):
        assert compare_results_no_decimal(
            fname,
            pd.read_csv(
                dir_anim_results / "matrix_aln_lengths_1.tab", index_col=0, sep="\t"
            ),
            0,
        )


def test_sim_errors(
    anim_nucmer_targets_filter_indir: Path,
    dir_anim_results: Path,
) -> None:
    """Check aniM sim errors."""
    for fname in (anim_nucmer_targets_filter_indir).glob("*.filter"):
        assert compare_results_no_decimal(
            fname,
            pd.read_csv(
                dir_anim_results / "matrix_sim_errors_1.tab", index_col=0, sep="\t"
            ),
            3,
        )


def test_delta_parsing(anim_nucmer_targets_filter_indir: Path) -> None:
    """Check parsing of test NUCmer .delta/.filter file."""
    assert method_anim.parse_delta(
        anim_nucmer_targets_filter_indir
        / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.filter",
    ) == (
        39169,
        39176,
        0.9962487643734,
        222,
    )


def test_aligned_bases_count(aligned_regions: dict) -> None:
    """Check only aligned bases in non-overlapping regions are counted."""
    assert method_anim.get_aligned_bases_count(aligned_regions) == 39176  # noqa: PLR2004


def test_get_genome_length(input_genomes_tiny: Path) -> None:
    """Check retrieved genome lengths match the expected values."""
    assert (
        method_anim.get_genome_length(input_genomes_tiny / "MGV-GENOME-0264574.fas")
        == 39253  # noqa: PLR2004
    )


def test_collect_results(
    expected_anim_results: method_anim.ComparisonResult,
    anim_nucmer_targets_filter_indir: Path,
    input_genomes_tiny: Path,
) -> None:
    """Check that ComparisionResult class is populated with correct values."""
    assert expected_anim_results == next(
        result
        for result in method_anim.collect_results_directory(
            anim_nucmer_targets_filter_indir, input_genomes_tiny
        )
        if result.qname == "MGV-GENOME-0264574" and result.rname == "MGV-GENOME-0266457"
    )
