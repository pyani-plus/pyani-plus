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
"""Test snakemake workflow for dnadiff.

These tests are intended to be run from the repository root using:

pytest -v or make test
"""

import shutil
from pathlib import Path

import pytest

from pyani_plus.snakemake import snakemake_scheduler
from pyani_plus.tools import (
    get_delta_filter,
    get_nucmer,
    get_show_coords,
    get_show_diff,
)


@pytest.fixture
def config_dnadiff_args(
    input_genomes_tiny: Path,
    snakemake_cores: int,
) -> dict:
    """Return configuration settings for testing snakemake dnadiff file.

    We take the output directories for the MUMmer delta/filter output and
    the small set of input genomes as arguments.
    """
    return {
        "nucmer": get_nucmer().exe_path,
        "delta_filter": get_delta_filter().exe_path,
        "show_coords": get_show_coords().exe_path,
        "show_diff": get_show_diff().exe_path,
        # "outdir": ... is dynamic
        "indir": str(input_genomes_tiny),
        "cores": snakemake_cores,
    }


def compare_files_with_skip(file1: Path, file2: Path, skip: int = 1) -> bool:
    """Compare two files, line by line, except for the first line.

    This function expects two text files as input and returns True if the content
    of the files is the same, and False if the two files differ.
    """
    with file1.open() as if1, file2.open() as if2:
        for line1, line2 in zip(
            if1.readlines()[skip:],
            if2.readlines()[skip:],
            strict=False,
        ):
            if line1 != line2:
                return False

    return True


def compare_show_diff_files(file1: Path, file2: Path) -> bool:
    """Compare two files.

    This function expects two text files as input and returns True if the content
    of the files is the same, and False if the two files differ.
    """
    with file1.open() as if1, file2.open() as if2:
        for line1, line2 in zip(if1.readlines(), if2.readlines(), strict=False):
            if line1 != line2:
                return False

    return True


def test_snakemake_rule_delta(
    dnadiff_nucmer_targets_delta: list[str],
    dnadiff_nucmer_targets_delta_indir: Path,
    dnadiff_nucmer_targets_delta_outdir: Path,
    config_dnadiff_args: dict,
    tmp_path: str,
) -> None:
    """Test nucmer delta snakemake wrapper.

    Checks that the delta rule in the dnadiff snakemake wrapper gives the
    expected output.

    If the output directory exists (i.e. the make clean_tests rule has not
    been run), the tests will automatically pass as snakemake will not
    attempt to re-run the rule. That would prevent us from seeing any
    introduced bugs, so we force re-running the rule by deleting the
    output directory before running the tests.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(dnadiff_nucmer_targets_delta_outdir, ignore_errors=True)

    config = config_dnadiff_args.copy()
    config["outdir"] = dnadiff_nucmer_targets_delta_outdir

    # Run snakemake wrapper
    runner = snakemake_scheduler.SnakemakeRunner("snakemake_dnadiff.smk")
    runner.run_workflow(dnadiff_nucmer_targets_delta, config, workdir=Path(tmp_path))

    # Check output against target fixtures
    for fname in dnadiff_nucmer_targets_delta:
        assert compare_files_with_skip(
            dnadiff_nucmer_targets_delta_indir / fname,
            dnadiff_nucmer_targets_delta_outdir / fname,
        )


def test_snakemake_rule_filter(
    dnadiff_nucmer_targets_filter: list[str],
    dnadiff_nucmer_targets_filter_indir: Path,
    dnadiff_nucmer_targets_filter_outdir: Path,
    config_dnadiff_args: dict,
    tmp_path: str,
) -> None:
    """Test nucmer filter snakemake wrapper.

    Checks that the filter rule in the dnadiff snakemake wrapper gives the
    expected output.

    If the output directory exists (i.e. the make clean_tests rule has not
    been run), the tests will automatically pass as snakemake will not
    attempt to re-run the rule. That would prevent us from seeing any
    introduced bugs, so we force re-running the rule by deleting the
    output directory before running the tests.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(dnadiff_nucmer_targets_filter_outdir, ignore_errors=True)

    config = config_dnadiff_args.copy()
    config["outdir"] = dnadiff_nucmer_targets_filter_outdir

    # Run snakemake wrapper
    runner = snakemake_scheduler.SnakemakeRunner("snakemake_dnadiff.smk")
    runner.run_workflow(dnadiff_nucmer_targets_filter, config, workdir=Path(tmp_path))

    # Check output against target fixtures
    for fname in dnadiff_nucmer_targets_filter:
        assert compare_files_with_skip(
            dnadiff_nucmer_targets_filter_indir / fname,
            dnadiff_nucmer_targets_filter_outdir / fname,
        )


def test_snakemake_rule_show_diff(
    dnadiff_targets_showdiff: list[str],
    dnadiff_targets_showdiff_indir: Path,
    dnadiff_targets_showdiff_outdir: Path,
    config_dnadiff_args: dict,
    tmp_path: str,
) -> None:
    """Test dnadiff show-diff snakemake wrapper.

    Checks that the show-diff rule in the dnadiff snakemake wrapper gives the
    expected output.

    If the output directory exists (i.e. the make clean_tests rule has not
    been run), the tests will automatically pass as snakemake will not
    attempt to re-run the rule. That would prevent us from seeing any
    introduced bugs, so we force re-running the rule by deleting the
    output directory before running the tests.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(dnadiff_targets_showdiff_outdir, ignore_errors=True)

    config = config_dnadiff_args.copy()
    config["outdir"] = dnadiff_targets_showdiff_outdir

    # Run snakemake wrapper
    runner = snakemake_scheduler.SnakemakeRunner("snakemake_dnadiff.smk")
    runner.run_workflow(dnadiff_targets_showdiff, config, workdir=Path(tmp_path))

    # Check output against target fixtures
    for fname in dnadiff_targets_showdiff:
        assert compare_files_with_skip(
            dnadiff_targets_showdiff_indir / fname,
            dnadiff_targets_showdiff_outdir / fname,
            skip=0,
        )


def test_snakemake_rule_show_coords(
    dnadiff_targets_showcoords: list[str],
    dnadiff_targets_showcoords_indir: Path,
    dnadiff_targets_showcoords_outdir: Path,
    config_dnadiff_args: dict,
    tmp_path: str,
) -> None:
    """Test dnadiff show-coords snakemake wrapper.

    Checks that the show-coords rule in the dnadiff snakemake wrapper gives the
    expected output.

    If the output directory exists (i.e. the make clean_tests rule has not
    been run), the tests will automatically pass as snakemake will not
    attempt to re-run the rule. That would prevent us from seeing any
    introduced bugs, so we force re-running the rule by deleting the
    output directory before running the tests.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(dnadiff_targets_showcoords_outdir, ignore_errors=True)

    config = config_dnadiff_args.copy()
    config["outdir"] = dnadiff_targets_showcoords_outdir

    # Run snakemake wrapper
    runner = snakemake_scheduler.SnakemakeRunner("snakemake_dnadiff.smk")
    runner.run_workflow(dnadiff_targets_showcoords, config, workdir=Path(tmp_path))

    # Check output against target fixtures
    for fname in dnadiff_targets_showcoords:
        assert compare_files_with_skip(
            dnadiff_targets_showcoords_indir / fname,
            dnadiff_targets_showcoords_outdir / fname,
            skip=0,
        )
