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
def dnadiff_nucmer_targets_filter_indir() -> Path:
    """Directory containing MUMmer filter snakemake reference files."""
    return FIXTUREPATH / "dnadiff" / "targets" / "filter"


@pytest.fixture
def dnadiff_nucmer_targets_delta_indir() -> Path:
    """Directory containing MUMmer delta snakemake reference files."""
    return FIXTUREPATH / "dnadiff" / "targets" / "delta"


@pytest.fixture
def dnadiff_targets_showdiff_indir() -> Path:
    """Directory containing snadiff show_diff snakemake reference files."""
    return FIXTUREPATH / "dnadiff" / "targets" / "show_diff"


@pytest.fixture
def dnadiff_targets_showcoords_indir() -> Path:
    """Directory containing snadiff show_diff snakemake reference files."""
    return FIXTUREPATH / "dnadiff" / "targets" / "show_coords"


@pytest.fixture
def dnadiff_targets_outdir(tmp_path: str) -> Path:
    """Output directory for MUMmer filter snakemake tests.

    This path indicates the location to which dnadiff should write
    its output files during dnadiff testing
    """
    return Path(tmp_path).resolve() / "dnadiff_targets_output"


@pytest.fixture
def dnadiff_nucmer_targets_filter(
    dnadiff_targets_outdir: Path,
) -> list[Path]:
    """Target filter files for dnadiff method tests.

    These are paths to the output files we want to generate using
    nucmer for dnadiff. We aim to ask MUMmer to generate a set of
    .filter files that could later be processed to obtain AlignedBases
    values
    """
    reference_paths = (FIXTUREPATH / "dnadiff" / "targets" / "filter").glob("*.filter")
    return [dnadiff_targets_outdir / _.name for _ in reference_paths]


@pytest.fixture
def dnadiff_nucmer_targets_delta(
    dnadiff_targets_outdir: Path,
) -> list[Path]:
    """Target delta files for dnadiff method tests.

    These are paths to the output files we want to generate using
    nucmer for dnadiff. We aim to generate a set of .delta files that
    could later be processed to obtain AlignedBases values
    """
    reference_paths = (FIXTUREPATH / "dnadiff" / "targets" / "delta").glob("*.delta")
    return [dnadiff_targets_outdir / _.name for _ in reference_paths]


@pytest.fixture
def dnadiff_targets_showdiff(dnadiff_targets_outdir: Path) -> list[Path]:
    """Target qdiff files for dnadiff method tests.

    These are paths to the output files we want to generate using
    show-diff for dnadiff method. We aim to generate a set of .qdiff files that
    could later be processed to obtain ANI values.
    """
    reference_paths = (FIXTUREPATH / "dnadiff" / "targets" / "show_diff").glob(
        "*.qdiff",
    )
    return [dnadiff_targets_outdir / _.name for _ in reference_paths]


@pytest.fixture
def dnadiff_targets_showcoords(dnadiff_targets_outdir: Path) -> list[Path]:
    """Target mcoords files for dnadiff method tests.

    These are paths to the output files we want to generate using
    show-coords for dnadiff method. We aim to generate a set of
    .mcoords files that could later be processed to obtain ANI values.
    """
    reference_paths = (FIXTUREPATH / "dnadiff" / "targets" / "show_coords").glob(
        "*.mcoords",
    )
    return [dnadiff_targets_outdir / _.name for _ in reference_paths]


@pytest.fixture
def input_genomes_tiny() -> Path:
    """Path to small set of two viral input genomes."""
    return FIXTUREPATH / "viral_example"


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
