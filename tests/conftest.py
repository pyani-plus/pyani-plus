# (c) The University of Strathclude 2024-present
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
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
# (c) The University of Strathclude 2024-present
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
"""Pytest configuration file."""

from pathlib import Path

import pytest

# Path to tests, contains tests and data subdirectories
# This conftest.py file should be found in the top directory of the tests
# module. The fixture data should be in a subdirectory named fixtures
TESTSPATH = Path(__file__).parents[0]
FIXTUREPATH = TESTSPATH / "fixtures"


@pytest.fixture
def anim_nucmer_targets_filter_outdir():
    """Output directory for MUMmer filter snakemake tests.

    This path indicates the location to which MUMmer should write
    its output files during ANIm testing
    """
    return TESTSPATH / "nucmer_filter_output"


@pytest.fixture
def anim_nucmer_targets_delta_outdir():
    """Output directory for MUMmer delta snakemake tests.

    This path indicates the location to which MUMmer should write
    its output files during ANIm testing
    """
    return TESTSPATH / "nucmer_delta_output"


@pytest.fixture
def anim_nucmer_targets_filter(anim_nucmer_targets_filter_outdir):
    """Target files for ANIm tests.

    These are paths to the output files we want to generate using
    nucmer for ANIm. We aim to ask MUMmer to generate a set of
    .filter files that could later be processed to obtain ANI values
    """
    reference_paths = (FIXTUREPATH / "anim" / "targets" / "filter").glob("*.filter")
    return [anim_nucmer_targets_filter_outdir / _.name for _ in reference_paths]


@pytest.fixture
def anim_nucmer_targets_delta(anim_nucmer_targets_delta_outdir):
    """Target files for ANIm tests.

    These are paths to the output files we want to generate using
    nucmer for ANIm. We aim to generate a set of .delta files that
    could later be processed to obtain ANI values
    """
    reference_paths = (FIXTUREPATH / "anim" / "targets" / "delta").glob("*.delta")
    return [anim_nucmer_targets_delta_outdir / _.name for _ in reference_paths]


@pytest.fixture
def anim_nucmer_expected_targets_filter():
    """Target files for ANIm tests.

    These are paths to the output files we want to generate using
    nucmer for ANIm. We aim to ask MUMmer to generate a set of
    .filter files that could later be processed to obtain ANI values
    """

    return FIXTUREPATH / "anim" / "targets" / "filter"


@pytest.fixture
def input_genomes_small():
    """Path to small set of input genomes."""
    return str(FIXTUREPATH / "sequences")
