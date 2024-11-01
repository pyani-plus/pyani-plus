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
"""Test snakemake workflow for sourmash.

These tests are intended to be run from the repository root using:

pytest -v or make test
"""

import shutil  # We need this for filesystem operations
from pathlib import Path

# Required to support pytest automated testing
import pytest

from pyani_plus.workflows import SnakemakeRunner


@pytest.fixture
def config_sourmash_args(
    input_genomes_tiny: Path,
    snakemake_cores: int,
) -> dict:
    """Return configuration settings for testing snakemake sourmash rules."""
    return {
        # "outdir": ... is dynamic
        "indir": input_genomes_tiny,
        "cores": snakemake_cores,
        "scaled": 1000,
        "kmer": 31,
    }


def compare_sourmash_files(file1: Path, file2: Path) -> bool:
    """Compare two files.

    This function expects two text files as input and returns True if the content
    of the files is the same, and False if the two files differ.
    """
    with file1.open() as if1, file2.open() as if2:
        for line1, line2 in zip(if1.readlines(), if2.readlines(), strict=False):
            if line1 != line2:
                return False

    return True


def test_snakemake_sketch_rule(
    sourmash_targets_signature_outdir: Path,
    sourmash_targets_sig: list[str],
    config_sourmash_args: dict,
    sourmash_targets_sig_indir: Path,
    tmp_path: str,
) -> None:
    """Test sourmash sketch snakemake wrapper.

    Checks that the sketch rule in the sourmash snakemake wrapper gives the
    expected output.

    If the output directory exists (i.e. the make clean_tests rule has not
    been run), the tests will automatically pass as snakemake will not
    attempt to re-run the rule. That would prevent us from seeing any
    introduced bugs, so we force re-running the rule by deleting the
    output directory before running the tests.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(sourmash_targets_signature_outdir, ignore_errors=True)

    config = config_sourmash_args.copy()
    config["outdir"] = sourmash_targets_signature_outdir

    # Run snakemake wrapper
    runner = SnakemakeRunner("snakemake_sourmash.smk")
    runner.run_workflow(sourmash_targets_sig, config, workdir=Path(tmp_path))

    # Check output against target fixtures
    # THIS TEST SHOULD FAIL, but passes instead?
    # Target/fixture files were generated with scaled=300
    # Here, I use scaled=1000?
    for fname in sourmash_targets_sig:
        assert compare_sourmash_files(
            sourmash_targets_sig_indir / fname,
            sourmash_targets_signature_outdir / fname,
        )


def test_snakemake_compare_rule(
    sourmash_targets_compare_outdir: Path,
    sourmash_compare_sig: list[str],
    config_sourmash_args: dict,
    sourmash_targets_compare_indir: Path,
    tmp_path: str,
) -> None:
    """Test sourmash compare snakemake wrapper.

    Checks that the compare rule in the sourmash snakemake wrapper gives the
    expected output.

    If the output directory exists (i.e. the make clean_tests rule has not
    been run), the tests will automatically pass as snakemake will not
    attempt to re-run the rule. That would prevent us from seeing any
    introduced bugs, so we force re-running the rule by deleting the
    output directory before running the tests.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(sourmash_targets_compare_outdir, ignore_errors=True)

    config = config_sourmash_args.copy()
    config["outdir"] = sourmash_targets_compare_outdir

    # Run snakemake wrapper
    runner = SnakemakeRunner("snakemake_sourmash.smk")
    runner.run_workflow(sourmash_compare_sig, config, workdir=Path(tmp_path))

    # Check output against target fixtures
    # THIS TEST SHOULD FAIL, but passes instead?
    # Target/fixture files were generated with scaled=300
    # Here, I use scaled=1000?
    for fname in sourmash_compare_sig:
        assert compare_sourmash_files(
            sourmash_targets_compare_indir / fname,
            sourmash_targets_compare_outdir / fname,
        )
