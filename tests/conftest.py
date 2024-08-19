# (c) The University of Strathclyde 2024-present
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
# (c) The University of Strathclyde 2024-present
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


@pytest.fixture()
def anim_nucmer_targets_filter_indir() -> Path:
    """Directory containing MUMmer filter snakemake reference files."""
    return FIXTUREPATH / "anim" / "targets" / "filter"


@pytest.fixture()
def anim_nucmer_targets_filter_outdir(tmp_path: str) -> Path:
    """Output directory for MUMmer filter snakemake tests.

    This path indicates the location to which MUMmer should write
    its output files during ANIm testing
    """
    return Path(tmp_path).resolve() / "nucmer_filter_output"


@pytest.fixture()
def anim_nucmer_targets_filter_slurm_outdir(tmp_path: str) -> Path:
    """Output directory for MUMmer filter snakemake tests.

    This path indicates the location to which MUMmer should write
    its output files during ANIm testing
    """
    return Path(tmp_path).resolve() / "nucmer_filter_slurm_output"


@pytest.fixture()
def anim_nucmer_targets_delta_indir() -> Path:
    """Directory containing MUMmer delta snakemake reference files."""
    return FIXTUREPATH / "anim" / "targets" / "delta"


@pytest.fixture()
def anim_nucmer_targets_delta_outdir(tmp_path: str) -> Path:
    """Output directory for MUMmer delta snakemake tests.

    This path indicates the location to which MUMmer should write
    its output files during ANIm testing
    """
    return Path(tmp_path).resolve() / "nucmer_delta_output"


@pytest.fixture()
def anim_nucmer_targets_delta_slurm_outdir(tmp_path: str) -> Path:
    """Output directory for MUMmer delta snakemake tests.

    This path indicates the location to which MUMmer should write
    its output files during ANIm testing
    """
    return Path(tmp_path).resolve() / "nucmer_delta_slurm_output"


@pytest.fixture()
def anim_nucmer_targets_filter(anim_nucmer_targets_filter_outdir: Path) -> list[str]:
    """Target files for ANIm tests.

    These are paths to the output files we want to generate using
    nucmer for ANIm. We aim to ask MUMmer to generate a set of
    .filter files that could later be processed to obtain ANI values
    """
    reference_paths = (FIXTUREPATH / "anim" / "targets" / "filter").glob("*.filter")
    return [anim_nucmer_targets_filter_outdir / _.name for _ in reference_paths]


@pytest.fixture()
def anim_nucmer_targets_filter_slurm(
    anim_nucmer_targets_filter_slurm_outdir: Path,
) -> list[str]:
    """Target files for ANIm tests.

    These are paths to the output files we want to generate using
    nucmer for ANIm. We aim to ask MUMmer to generate a set of
    .filter files that could later be processed to obtain ANI values
    """
    reference_paths = (FIXTUREPATH / "anim" / "targets" / "filter").glob("*.filter")
    return [anim_nucmer_targets_filter_slurm_outdir / _.name for _ in reference_paths]


@pytest.fixture()
def anim_nucmer_targets_delta(anim_nucmer_targets_delta_outdir: Path) -> list[str]:
    """Target files for ANIm tests.

    These are paths to the output files we want to generate using
    nucmer for ANIm. We aim to generate a set of .delta files that
    could later be processed to obtain ANI values
    """
    reference_paths = (FIXTUREPATH / "anim" / "targets" / "delta").glob("*.delta")
    return [anim_nucmer_targets_delta_outdir / _.name for _ in reference_paths]


@pytest.fixture()
def anim_nucmer_targets_delta_slurm(
    anim_nucmer_targets_delta_slurm_outdir: Path,
) -> list[str]:
    """Target files for ANIm tests.

    These are paths to the output files we want to generate using
    nucmer for ANIm. We aim to generate a set of .delta files that
    could later be processed to obtain ANI values
    """
    reference_paths = (FIXTUREPATH / "anim" / "targets" / "delta").glob("*.delta")
    return [anim_nucmer_targets_delta_slurm_outdir / _.name for _ in reference_paths]


@pytest.fixture()
def fastani_targets_outdir(tmp_path: str) -> Path:
    """Output directory for fastani snakemake tests.

    This path indicates the location to which fastANI should write
    its output files during fastani testing
    """
    return Path(tmp_path).resolve() / "fastani_output"


@pytest.fixture()
def fastani_targets_slurm_outdir(tmp_path: str) -> Path:
    """Output directory for fastani snakemake tests.

    This path indicates the location to which fastANI should write
    its output files during fastani testing
    """
    return Path(tmp_path).resolve() / "fastani_slurm_output"


@pytest.fixture()
def fastani_targets(fastani_targets_outdir: Path) -> list[str]:
    """Target files for ANIm tests.

    These are paths to the output files we want to generate using
    nucmer for ANIm. We aim to ask MUMmer to generate a set of
    .filter files that could later be processed to obtain ANI values
    """
    reference_paths = (FIXTUREPATH / "fastani" / "targets").glob("*.fastani")
    return [fastani_targets_outdir / _.name for _ in reference_paths]


@pytest.fixture()
def fastani_targets_slurm(fastani_targets_slurm_outdir: Path) -> list[str]:
    """Target files for ANIm tests.

    These are paths to the output files we want to generate using
    nucmer for ANIm. We aim to ask MUMmer to generate a set of
    .filter files that could later be processed to obtain ANI values
    """
    reference_paths = (FIXTUREPATH / "fastani" / "targets").glob("*.fastani")
    return [fastani_targets_slurm_outdir / _.name for _ in reference_paths]


@pytest.fixture()
def fastani_targets_indir() -> Path:
    """Target files for fastANI tests.

    These are paths to the output files we generated using
    fastani. We aim to use these to test if the fastani
    files generated by snakemake workflow are the same.
    """
    return FIXTUREPATH / "fastani" / "targets"


@pytest.fixture()
def dnadiff_nucmer_targets_filter_indir() -> Path:
    """Directory containing MUMmer filter snakemake reference files."""
    return FIXTUREPATH / "dnadiff" / "targets" / "filter"


@pytest.fixture()
def dnadiff_nucmer_targets_filter_outdir(tmp_path: str) -> Path:
    """Output directory for MUMmer filter snakemake tests.

    This path indicates the location to which MUMmer should write
    its output files during dnadiff testing
    """
    return Path(tmp_path).resolve() / "dnadiff_nucmer_filter_output"


@pytest.fixture()
def dnadiff_nucmer_targets_filter_slurm_outdir(tmp_path: str) -> Path:
    """Output directory for MUMmer filter snakemake tests.

    This path indicates the location to which MUMmer should write
    its output files during dnadiff testing
    """
    return Path(tmp_path).resolve() / "dnadiff_nucmer_filter_slurm_output"


@pytest.fixture()
def dnadiff_nucmer_targets_delta_indir() -> Path:
    """Directory containing MUMmer delta snakemake reference files."""
    return FIXTUREPATH / "dnadiff" / "targets" / "delta"


@pytest.fixture()
def dnadiff_nucmer_targets_delta_outdir(tmp_path: str) -> Path:
    """Output directory for MUMmer filter snakemake tests.

    This path indicates the location to which MUMmer should write
    its output files during dnadiff testing
    """
    return Path(tmp_path).resolve() / "dnadiff_nucmer_delta_output"


@pytest.fixture()
def dnadiff_nucmer_targets_delta_slurm_outdir(tmp_path: str) -> Path:
    """Output directory for MUMmer filter snakemake tests.

    This path indicates the location to which MUMmer should write
    its output files during dnadiff testing
    """
    return Path(tmp_path).resolve() / "dnadiff_nucmer_delta_slurm_output"


@pytest.fixture()
def dnadiff_targets_showdiff_indir() -> Path:
    """Directory containing snadiff show_diff snakemake reference files."""
    return FIXTUREPATH / "dnadiff" / "targets" / "show_diff"


@pytest.fixture()
def dnadiff_targets_showdiff_outdir(tmp_path: str) -> Path:
    """Output directory for MUMmer filter snakemake tests.

    This path indicates the location to which dnadiff should write
    its output files during dnadiff testing
    """
    return Path(tmp_path).resolve() / "dnadiff_show_diff_output"


@pytest.fixture()
def dnadiff_targets_showcoords_indir() -> Path:
    """Directory containing snadiff show_diff snakemake reference files."""
    return FIXTUREPATH / "dnadiff" / "targets" / "show_coords"


@pytest.fixture()
def dnadiff_targets_showcoords_outdir(tmp_path: str) -> Path:
    """Output directory for MUMmer filter snakemake tests.

    This path indicates the location to which dnadiff should write
    its output files during dnadiff testing
    """
    return Path(tmp_path).resolve() / "dnadiff_show_coords_output"


@pytest.fixture()
def dnadiff_targets_showdiff_slurm_outdir(tmp_path: str) -> Path:
    """Output directory for MUMmer filter snakemake tests.

    This path indicates the location to which dnadiff should write
    its output files during dnadiff testing
    """
    return Path(tmp_path).resolve() / "dnadiff_showdiff_slurm_output"


@pytest.fixture()
def dnadiff_nucmer_targets_filter(
    dnadiff_nucmer_targets_filter_outdir: Path,
) -> list[str]:
    """Target files for dnadiff tests.

    These are paths to the output files we want to generate using
    nucmer for dnadiff. We aim to ask MUMmer to generate a set of
    .filter files that could later be processed to obtain AlignedBases
    values
    """
    reference_paths = (FIXTUREPATH / "dnadiff" / "targets" / "filter").glob("*.filter")
    return [dnadiff_nucmer_targets_filter_outdir / _.name for _ in reference_paths]


@pytest.fixture()
def dnadiff_nucmer_targets_filter_slurm(
    dnadiff_nucmer_targets_filter_slurm_outdir: Path,
) -> list[str]:
    """Target files for dnadiff tests.

    These are paths to the output files we want to generate using
    nucmer for dnadiff. We aim to ask MUMmer to generate a set of
    .filter files that could later be processed to obtain AlignedBases
    values
    """
    reference_paths = (FIXTUREPATH / "dnadiff" / "targets" / "filter").glob("*.filter")
    return [
        dnadiff_nucmer_targets_filter_slurm_outdir / _.name for _ in reference_paths
    ]


@pytest.fixture()
def dnadiff_nucmer_targets_delta(
    dnadiff_nucmer_targets_delta_outdir: Path,
) -> list[str]:
    """Target files for dnadiff tests.

    These are paths to the output files we want to generate using
    nucmer for dnadiff. We aim to generate a set of .delta files that
    could later be processed to obtain AlignedBases values
    """
    reference_paths = (FIXTUREPATH / "dnadiff" / "targets" / "delta").glob("*.delta")
    return [dnadiff_nucmer_targets_delta_outdir / _.name for _ in reference_paths]


@pytest.fixture()
def dnadiff_nucmer_targets_delta_slurm(
    dnadiff_nucmer_targets_delta_slurm_outdir: Path,
) -> list[str]:
    """Target files for dnadiff tests.

    These are paths to the output files we want to generate using
    nucmer for dnadiff. We aim to generate a set of .delta files that
    could later be processed to obtain AlignedBases values
    """
    reference_paths = (FIXTUREPATH / "dnadiff" / "targets" / "delta").glob("*.delta")
    return [dnadiff_nucmer_targets_delta_slurm_outdir / _.name for _ in reference_paths]


@pytest.fixture()
def dnadiff_targets_showdiff(dnadiff_targets_showdiff_outdir: Path) -> list[str]:
    """Target files for ANIm tests.

    These are paths to the output files we want to generate using
    nucmer for ANIm. We aim to generate a set of .delta files that
    could later be processed to obtain ANI values
    """
    reference_paths = (FIXTUREPATH / "dnadiff" / "targets" / "show_diff").glob(
        "*.rdiff",
    )
    return [dnadiff_targets_showdiff_outdir / _.name for _ in reference_paths]


@pytest.fixture()
def dnadiff_targets_showcoords(dnadiff_targets_showcoords_outdir: Path) -> list[Path]:
    """Target files for ANIm tests.

    These are paths to the output files we want to generate using
    nucmer for ANIm. We aim to generate a set of .delta files that
    could later be processed to obtain ANI values
    """
    reference_paths = (FIXTUREPATH / "dnadiff" / "targets" / "show_coords").glob(
        "*.mcoords",
    )
    return [dnadiff_targets_showcoords_outdir / _.name for _ in reference_paths]


@pytest.fixture()
def dnadiff_targets_showdiff_slurm(
    dnadiff_targets_showdiff_slurm_outdir: Path,
) -> list[str]:
    """Target files for ANIm tests.

    These are paths to the output files we want to generate using
    nucmer for ANIm. We aim to generate a set of .delta files that
    could later be processed to obtain ANI values
    """
    reference_paths = (FIXTUREPATH / "dnadiff" / "targets" / "show_diff").glob(
        "*.rdiff",
    )
    return [dnadiff_targets_showdiff_slurm_outdir / _.name for _ in reference_paths]


@pytest.fixture()
def input_genomes_small() -> Path:
    """Path to small set of input genomes."""
    return FIXTUREPATH / "sequences"
