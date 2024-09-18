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
"""Test snakemake workflow for indexing genomes (MD5 checksums etc).

These tests are intended to be run from the repository root using:

    $ make test

Or,

    $ pytest tests/snakemake/test_snakemake_index_workflow.py
"""

import shutil
from pathlib import Path

import pytest

from pyani_plus.snakemake import snakemake_scheduler


def compare_checksum_files(file1: Path, file2: Path) -> bool:
    """Compare contents of files (here MD5 checksum files)."""
    with file1.open("rb") as handle:
        data1 = handle.read()
    with file2.open("rb") as handle:
        data2 = handle.read()
    return data1 == data2


@pytest.fixture
def config_fastani_args(
    index_targets_outdir: Path,
    input_genomes_small: Path,
    snakemake_cores: int,
) -> dict:
    """Return configuration settings for testing snakemake fastANI rule."""
    return {
        "md5sum": "md5sum",
        "outdir": index_targets_outdir,
        "indir": str(input_genomes_small),
        "cores": snakemake_cores,
    }


def test_snakemake_rule_index(
    index_targets: list[str],
    index_targets_indir: Path,
    index_targets_outdir: Path,
    config_fastani_args: dict,
    tmp_path: str,
) -> None:
    """Test index snakemake wrapper.

    Checks that the index snakemake wrapper gives the expected output.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(index_targets_outdir, ignore_errors=True)

    # Run snakemake wrapper
    runner = snakemake_scheduler.SnakemakeRunner("snakemake_index.smk")

    runner.run_workflow(index_targets, config_fastani_args, workdir=Path(tmp_path))

    # Check output against target fixtures
    for fname in index_targets:
        assert compare_checksum_files(
            index_targets_indir / fname,
            index_targets_outdir / fname,
        )
