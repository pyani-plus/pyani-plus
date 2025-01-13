# The MIT License
#
# Copyright (c) 2024-2025 University of Strathclyde
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
# Resolve this to an absolute path so that the working directory can be changed
TESTSPATH = Path(__file__).parents[0].resolve()
FIXTUREPATH = TESTSPATH / "fixtures"


@pytest.fixture
def anib_targets_outdir(tmp_path: str) -> Path:
    """Output directory for ANIb snakemake tests."""
    return Path(tmp_path).resolve() / "anib_output"


@pytest.fixture
def anim_targets_outdir(tmp_path: str) -> Path:
    """Output directory for ANIm snakemake tests.

    This path indicates the location to which MUMmer should write
    its output files during ANIm testing
    """
    return Path(tmp_path).resolve() / "nucmer_filter_output"


@pytest.fixture
def fastani_targets_outdir(tmp_path: str) -> Path:
    """Output directory for fastani snakemake tests.

    This path indicates the location to which fastANI should write
    its output files during fastani testing
    """
    return Path(tmp_path).resolve() / "fastani_output"


@pytest.fixture
def dnadiff_targets_outdir(tmp_path: str) -> Path:
    """Output directory for MUMmer filter snakemake tests.

    This path indicates the location to which dnadiff should write
    its output files during dnadiff testing
    """
    return Path(tmp_path).resolve() / "dnadiff_targets_output"


@pytest.fixture
def sourmash_targets_signature_outdir(tmp_path: str) -> Path:
    """Output directory for sourmash sketch_dna snakemake tests.

    This path indicates the location to which sourmash should write
    its output files during sourmash testing
    """
    return Path(tmp_path).resolve() / "sourmash_sketch_output"


@pytest.fixture
def sourmash_targets_compare_outdir(tmp_path: str) -> Path:
    """Output directory for sourmash sketch_dna snakemake tests.

    This path indicates the location to which sourmash should write
    its output files during sourmash testing
    """
    return Path(tmp_path).resolve() / "sourmash_compare_output"


@pytest.fixture(scope="session")
def input_genomes_tiny() -> Path:
    """Path to small set of three viral input genomes."""
    return FIXTUREPATH / "viral_example"


@pytest.fixture(scope="session")
def input_genomes_bad_alignments() -> Path:
    """Path to small set of two bad alignments input genomes."""
    return FIXTUREPATH / "bad_alignments"


@pytest.fixture(scope="session")
def input_gzip_bacteria() -> Path:
    """Path to small set of four gzipped bacterial input genomes."""
    return FIXTUREPATH / "bacterial_example"
