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

import filecmp
import shutil  # We need this for filesystem operations
from pathlib import Path

# Required to support pytest automated testing
import pytest

# We're testing the workflow, as called through the pyani_plus
# wrapper code, so import only that module
from pyani_plus.snakemake import snakemake_scheduler
from pyani_plus.tools import get_blastn, get_makeblastdb


@pytest.fixture
def config_anib_args(
    anib_targets_outdir: Path,
    input_genomes_tiny: Path,
    snakemake_cores: int,
) -> dict:
    """Return configuration settings for testing snakemake filter rule.

    We take the output directories for the MUMmer filter output and the
    small set of input genomes as arguments.
    """
    return {
        "blastn": get_blastn().exe_path,
        "makeblastdb": get_makeblastdb().exe_path,
        "outdir": anib_targets_outdir,
        "indir": input_genomes_tiny,
        "cores": snakemake_cores,
        "fragLen": "1020",
    }


def test_snakemake_rule_fragments(
    anib_targets_fragments: list[Path],
    anib_targets_outdir: Path,
    config_anib_args: dict,
    anib_fragments: Path,
    tmp_path: str,
) -> None:
    """Test overall ANIb snakemake wrapper."""
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(anib_targets_outdir, ignore_errors=True)

    # Run snakemake wrapper
    runner = snakemake_scheduler.SnakemakeRunner("snakemake_anib.smk")
    runner.run_workflow(
        anib_targets_fragments, config_anib_args, workdir=Path(tmp_path)
    )

    # Check output against target fixtures
    for fname in anib_targets_fragments:
        assert fname.name.endswith("-fragments.fna"), fname
        assert Path(fname).is_file()
        assert filecmp.cmp(fname, anib_fragments / fname.name)


def compare_blast_json(file_a: Path, file_b: Path) -> bool:
    """Compare two BLAST+ .njs JSON files, ignoring the date-stamp."""
    with file_a.open() as handle_a, file_b.open() as handle_b:
        for a, b in zip(handle_a, handle_b, strict=True):
            assert (
                a == b
                or ("last-updated" in a and "last-updated" in b)
                or ("bytes-total" in a and "bytes-total" in b)
                or ("bytes-to-cache" in a and "bytes-to-cache" in b)
            ), f"{a!r} != {b!r}"
    return True


def test_snakemake_rule_blastdb(
    anib_targets_blastdb: list[Path],
    anib_targets_outdir: Path,
    config_anib_args: dict,
    anib_blastdb: Path,
    tmp_path: str,
) -> None:
    """Test overall ANIb snakemake wrapper."""
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(anib_targets_outdir, ignore_errors=True)

    # Run snakemake wrapper
    runner = snakemake_scheduler.SnakemakeRunner("snakemake_anib.smk")
    runner.run_workflow(anib_targets_blastdb, config_anib_args, workdir=Path(tmp_path))

    # Check output against target fixtures
    for fname in anib_targets_blastdb:
        assert fname.suffix == ".njs", fname
        assert Path(fname).is_file()
        compare_blast_json(fname, anib_blastdb / fname.name)


def test_snakemake_rule_blastn(
    anib_targets_blastn: list[Path],
    anib_targets_outdir: Path,
    config_anib_args: dict,
    anib_blastn: Path,
    tmp_path: str,
) -> None:
    """Test blastn (overall) ANIb snakemake wrapper."""
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(anib_targets_outdir, ignore_errors=True)

    # Run snakemake wrapper
    runner = snakemake_scheduler.SnakemakeRunner("snakemake_anib.smk")
    runner.run_workflow(anib_targets_blastn, config_anib_args, workdir=Path(tmp_path))

    # Check output against target fixtures
    for fname in anib_targets_blastn:
        assert fname.suffix == ".tsv", fname
        assert Path(fname).is_file()
        assert (anib_blastn / fname.name).is_file()
        assert filecmp.cmp(fname, anib_blastn / fname.name)
