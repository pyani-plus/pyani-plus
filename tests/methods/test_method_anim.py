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
        "MGV_MGV-GENOME-0266457": [(1, 37636), (17626, 39176)],
    }


@pytest.fixture
def anim_results(dir_anim_results: Path) -> pd.DataFrame:
    """Return expected results for anim."""
    return pd.read_csv(
        dir_anim_results / "matrix_identity_1.tab", index_col=0, sep="\t"
    )


def compare_anim_results(filterfile: Path, datadir: pd.DataFrame) -> bool:
    """Compare delta-filter files and expected anim results.

    This function expects delta-filter file and dataframe with expected
    anim results as input and returns True if the anim values is the same,
    and False if the two files differ.
    """
    reference, query = filterfile.stem.split("_vs_")

    return (
        round(method_anim.parse_delta(filterfile)[2], 6)
        == datadir.loc[reference, query]
    )


def test_anim_average_identity(
    anim_nucmer_targets_filter_indir: Path,
    anim_results: pd.DataFrame,
) -> None:
    """Check aniM average identity."""
    for fname in (
        anim_nucmer_targets_filter_indir / "anim" / "targets" / "filter"
    ).glob("*.filter"):
        assert compare_anim_results(fname, anim_results)


def test_anim_parsing(anim_nucmer_targets_filter_indir: Path) -> None:
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
    """Check aligned bases for aniM method."""
    assert method_anim.get_aligned_bases_count(aligned_regions) == 39176  # noqa: PLR2004


def test_get_genome_length(input_genomes_tiny: Path) -> None:
    """Check if genome lengths are correct."""
    assert (
        method_anim.get_genome_length(input_genomes_tiny / "MGV-GENOME-0264574.fna")
        == 39253  # noqa: PLR2004
    )


# TODO @kiepczi: Add a test for method_anim.py::collect_results_directory()
# https://github.com/pyani-plus/pyani-plus/issues/4
