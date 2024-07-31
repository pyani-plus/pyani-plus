# (c) University of Strathclyde 2019-
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2024-present University of Strathclyde
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
"""Test snakemake workflow for fastANI.

These tests are intended to be run from the repository root using:

make test
"""

import shutil
from pathlib import Path

import pytest

from pyani_plus.snakemake import fastani


@pytest.fixture()
def config_fastani_args(
    fastani_targets_outdir: Path,
    input_genomes_small: Path,
) -> dict:
    """Return configuration settings for testing snakemake fastANI rule."""
    return {
        "outdir": fastani_targets_outdir,
        "indir": str(input_genomes_small),
        "cores": 8,
        "fragLen": 3000,
        "kmerSize": 16,
        "minFrac": 0.2,
    }


def compare_fastani_files(file1: Path, file2: Path) -> bool:
    """Compare two fastANI files.

    This function expects two text files as input and returns True if the content
    of the files is the same, and False if the two files differ.

    As the Path to both Query and Reference might be diffrent,
    we will only consider file name.
    """
    with file1.open() as if1, file2.open() as if2:
        try:
            for line1, line2 in zip(if1.readlines(), if2.readlines(), strict=True):
                if line1.split("\t")[2:] != line2.split("\t")[2:]:
                    return False
            return True
        except ValueError:
            return False


def test_snakemake_rule_fastani(
    fastani_targets: list[str],
    fastani_targets_outdir: Path,
    config_fastani_args: dict,
    fastani_targets_indir: Path,
) -> None:
    """Test fastANI snakemake wrapper.

    Checks that the fastANI rule in the fastANI snakemake wrapper gives the
    expected output.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(fastani_targets_outdir, ignore_errors=True)

    # Run snakemake wrapper
    fastani.run_workflow(fastani_targets, config_fastani_args)

    # Check output against target fixtures
    for fname in fastani_targets:
        assert compare_fastani_files(  # noqa: S101
            fastani_targets_indir / fname,
            fastani_targets_outdir / fname,
        )
