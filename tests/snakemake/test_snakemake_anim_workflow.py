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
"""Test snakemake workflow for ANIm.

These tests are intended to be run from the repository root using:

pytest -v or make test
"""

import shutil  # We need this for filesystem operations
from pathlib import Path

# Required to support pytest automated testing
import pytest

from pyani_plus.private_cli import log_configuration, log_genome, log_run
from pyani_plus.snakemake import snakemake_scheduler
from pyani_plus.tools import get_delta_filter, get_nucmer

from . import compare_matrices


@pytest.fixture
def config_anim_args(
    input_genomes_tiny: Path,
    snakemake_cores: int,
    tmp_path: str,
) -> dict:
    """Return configuration settings for testing snakemake filter rule.

    We take the output directories for the MUMmer filter output and the
    small set of input genomes as arguments.
    """
    return {
        "db": Path(tmp_path) / "db.sqlite",
        "nucmer": get_nucmer().exe_path,
        "delta_filter": get_delta_filter().exe_path,
        # "outdir": ... is dynamic
        "indir": input_genomes_tiny,
        "cores": snakemake_cores,
    }


def compare_files_with_skip(file1: Path, file2: Path, skip: int = 1) -> bool:
    """Compare two files, line by line, skipping an initial count of lines.

    - skip: int, number of initial lines to skip in each file

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


def test_snakemake_rule_filter(  # noqa: PLR0913
    anim_nucmer_targets_filter: list[str],
    anim_nucmer_targets_filter_indir: Path,
    anim_nucmer_targets_filter_outdir: Path,
    dir_anim_results: Path,
    config_anim_args: dict,
    tmp_path: str,
) -> None:
    """Test nucmer filter snakemake wrapper.

    Checks that the filter rule in the ANIm snakemake wrapper gives the
    expected output.

    If the output directory exists (i.e. the make clean_tests rule has not
    been run), the tests will automatically pass as snakemake will not
    attempt to re-run the rule. That would prevent us from seeing any
    introduced bugs, so we force re-running the rule by deleting the
    output directory before running the tests.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(anim_nucmer_targets_filter_outdir, ignore_errors=True)

    # Assuming this will match but worker nodes might have a different version
    nucmer_tool = get_nucmer()

    # Setup minimal test DB
    db = config_anim_args["db"]
    assert not db.is_file()
    log_configuration(
        database=db,
        method="ANIm",
        program=nucmer_tool.exe_path.stem,
        version=nucmer_tool.version,
        fragsize=None,
        maxmatch=None,
        kmersize=None,
        minmatch=None,
        create_db=True,
    )
    # Record the FASTA files in the genomes table _before_ call snakemake
    log_genome(
        database=db,
        fasta=list(
            snakemake_scheduler.check_input_stems(config_anim_args["indir"]).values()
        ),
    )
    assert db.is_file()

    config = config_anim_args.copy()
    config["outdir"] = anim_nucmer_targets_filter_outdir

    # Run snakemake wrapper
    runner = snakemake_scheduler.SnakemakeRunner("snakemake_anim.smk")
    runner.run_workflow(anim_nucmer_targets_filter, config, workdir=Path(tmp_path))

    # Check output against target fixtures
    for fname in anim_nucmer_targets_filter:
        assert compare_files_with_skip(
            anim_nucmer_targets_filter_indir / fname,
            anim_nucmer_targets_filter_outdir / fname,
        )

    log_run(
        fasta=config_anim_args["indir"].glob("*.f*"),
        database=db,
        status="Complete",
        name="Test case",
        cmdline="pyani-plus anib --database ... blah blah blah",
        method="ANIm",
        program=nucmer_tool.exe_path.stem,
        version=nucmer_tool.version,
        fragsize=None,
        maxmatch=None,
        kmersize=None,
        minmatch=None,
        create_db=False,
    )
    compare_matrices(db, dir_anim_results)


def test_snakemake_rule_delta(
    anim_nucmer_targets_delta: list[str],
    anim_nucmer_targets_delta_indir: Path,
    anim_nucmer_targets_delta_outdir: Path,
    config_anim_args: dict,
    tmp_path: str,
) -> None:
    """Test nucmer delta snakemake wrapper.

    Checks that the delta rule in the ANIm snakemake wrapper gives the
    expected output.

    If the output directory exists (i.e. the make clean_tests rule has not
    been run), the tests will automatically pass as snakemake will not
    attempt to re-run the rule. That would prevent us from seeing any
    introduced bugs, so we force re-running the rule by deleting the
    output directory before running the tests.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(anim_nucmer_targets_delta_outdir, ignore_errors=True)

    config = config_anim_args.copy()
    config["outdir"] = anim_nucmer_targets_delta_outdir

    # Run snakemake wrapper
    runner = snakemake_scheduler.SnakemakeRunner("snakemake_anim.smk")
    runner.run_workflow(anim_nucmer_targets_delta, config, workdir=Path(tmp_path))

    # Check output against target fixtures
    for fname in anim_nucmer_targets_delta:
        assert compare_files_with_skip(
            anim_nucmer_targets_delta_indir / fname,
            anim_nucmer_targets_delta_outdir / fname,
        )
