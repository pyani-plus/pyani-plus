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

from pyani_plus.private_cli import log_run
from pyani_plus.tools import get_delta_filter, get_nucmer
from pyani_plus.workflows import (
    ToolExecutor,
    run_snakemake_with_progress_bar,
)

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
        "run_id": 1,  # by construction
        "nucmer": get_nucmer().exe_path,
        "delta_filter": get_delta_filter().exe_path,
        # "outdir": ... is dynamic
        "indir": input_genomes_tiny,
        "mode": "mum",
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


def test_rule_ANIm(  # noqa: N802, PLR0913
    input_genomes_tiny: Path,
    anim_nucmer_targets_filter: list[Path],
    anim_nucmer_targets_delta: list[Path],
    anim_nucmer_targets_filter_indir: Path,
    anim_targets_outdir: Path,
    anim_nucmer_targets_delta_indir: Path,
    config_anim_args: dict,
    tmp_path: str,
) -> None:
    """Test rule.

    Checks that the ANIm rule in the ANIm snakemake wrapper gives the
    expected output.

    If the output directory exists (i.e. the make clean_tests rule has not
    been run), the tests will automatically pass as snakemake will not
    attempt to re-run the rule. That would prevent us from seeing any
    introduced bugs, so we force re-running the rule by deleting the
    output directory before running the tests.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(anim_targets_outdir, ignore_errors=True)

    # Assuming this will match but worker nodes might have a different version
    nucmer_tool = get_nucmer()

    # Setup minimal test DB
    db = config_anim_args["db"]
    assert not db.is_file()
    log_run(
        fasta=config_anim_args["indir"],  # i.e. input_genomes_tiny
        database=db,
        status="Testing",
        name="Test case",
        cmdline="pyani-plus anib --database ... blah blah blah",
        method="ANIm",
        program=nucmer_tool.exe_path.stem,
        version=nucmer_tool.version,
        mode=config_anim_args["mode"],
        create_db=True,
    )
    assert db.is_file()

    config = config_anim_args.copy()
    config["outdir"] = anim_targets_outdir

    # Run snakemake wrapper
    run_snakemake_with_progress_bar(
        executor=ToolExecutor.local,
        workflow_name="snakemake_anim.smk",
        database=db,
        run_id=0,  # only needed for progress bar
        targets=anim_nucmer_targets_filter,
        params=config,
        working_directory=Path(tmp_path),
    )

    # Check delta-filter output against target fixtures
    for fname in anim_nucmer_targets_filter:
        assert compare_files_with_skip(
            anim_nucmer_targets_filter_indir / fname.name,
            fname,
        )

    # Check nucmer output against target fixtures
    for fname in anim_nucmer_targets_delta:
        assert compare_files_with_skip(
            anim_nucmer_targets_delta_indir / fname.name,
            fname,
        )

    compare_matrices(db, input_genomes_tiny / "matrices", "ANIm")
