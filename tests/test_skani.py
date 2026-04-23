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
"""Test methods for calculating ANI using skani.

These tests are intended to be run from the repository root using:

make test
"""

from pathlib import Path

from pyani_plus.methods import skani


def test_skani_parsing(input_genomes_tiny: Path) -> None:
    """Check parsing of test skani output files."""
    assert skani.parse_skani(
        input_genomes_tiny
        / "intermediates"
        / "skani"
        / "689d3fd6881db36b5e08329cf23cecdd_vs_689d3fd6881db36b5e08329cf23cecdd.skani"
    ) == (
        "../fixtures/viral_example/MGV-GENOME-0264574.fas",
        "../fixtures/viral_example/MGV-GENOME-0264574.fas",
        100.00,
        99.65,
        99.65,
        "MGV_MGV-GENOME-0264574",
        "MGV_MGV-GENOME-0264574",
    )

    assert skani.parse_skani(
        input_genomes_tiny
        / "intermediates"
        / "skani"
        / "78975d5144a1cd12e98898d573cf6536_vs_5584c7029328dc48d33f95f0a78f7e57.skani"
    ) == (
        "../fixtures/viral_example/OP073605.fasta",
        "../fixtures/viral_example/MGV-GENOME-0266457.fna",
        99.64,
        68.28,
        99.66,
        "OP073605.1 MAG: Bacteriophage sp. isolate 0984_12761, complete genome",
        "MGV_MGV-GENOME-0266457",
    )
