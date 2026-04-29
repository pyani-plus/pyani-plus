# The MIT License
#
# Copyright (c) 2024-2026 University of Strathclyde
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
"""Test methods for calculating ANI using lz-ani.

These tests are intended to be run from the repository root using:

make test
"""

from pathlib import Path

import pytest

from pyani_plus.methods import lzani


@pytest.fixture
def lzani_bad_headers() -> Path:
    """Path to lz-ani output file with bad headers."""
    return Path(__file__).parent / "fixtures" / "lzani" / "bad_headers.tsv"


def test_lzani_parsing(input_genomes_tiny: Path) -> None:
    """Check parsing of test lz-ani output files."""
    assert lzani.parse_lzani(
        input_genomes_tiny
        / "intermediates"
        / "lzani"
        / "689d3fd6881db36b5e08329cf23cecdd_vs_5584c7029328dc48d33f95f0a78f7e57.tsv"
    ) == (
        ("MGV-GENOME-0264574.fas", "OP073605.fasta", 0.999413, 0.998268),
        ("OP073605.fasta", "MGV-GENOME-0264574.fas", 0.987119, 0.701227),
    )


def test_lzani_parsing_bad_headers(lzani_bad_headers: Path) -> None:
    """Check parsing of lz-ani output files with bad headers."""
    with pytest.raises(SystemExit, match="Unexpected lz-ani output file format"):
        lzani.parse_lzani(lzani_bad_headers)
