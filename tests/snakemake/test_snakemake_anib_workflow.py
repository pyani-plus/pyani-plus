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

from pyani_plus.private_cli import log_configuration, log_genome, log_run
from pyani_plus.tools import get_blastn, get_makeblastdb
from pyani_plus.workflows import SnakemakeRunner, check_input_stems

from . import compare_matrices


@pytest.fixture
def config_anib_args(
    anib_targets_outdir: Path,
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
        "blastn": get_blastn().exe_path,
        "makeblastdb": get_makeblastdb().exe_path,
        "outdir": anib_targets_outdir,
        "indir": input_genomes_tiny,
        "cores": snakemake_cores,
        "fragsize": "1020",
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
    runner = SnakemakeRunner("snakemake_anib.smk")
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
    runner = SnakemakeRunner("snakemake_anib.smk")
    runner.run_workflow(anib_targets_blastdb, config_anib_args, workdir=Path(tmp_path))

    # Check output against target fixtures
    for fname in anib_targets_blastdb:
        assert fname.suffix == ".njs", fname
        assert Path(fname).is_file()
        compare_blast_json(fname, anib_blastdb / fname.name)


def test_snakemake_rule_blastn(  # noqa: PLR0913
    anib_targets_blastn: list[Path],
    anib_targets_outdir: Path,
    config_anib_args: dict,
    anib_blastn: Path,
    dir_anib_results: Path,
    tmp_path: str,
) -> None:
    """Test blastn (overall) ANIb snakemake wrapper."""
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(anib_targets_outdir, ignore_errors=True)

    # Assuming this will match but worker nodes might have a different version
    blastn_tool = get_blastn()

    # Setup minimal test DB
    db = config_anib_args["db"]
    assert not db.is_file()
    log_configuration(
        database=db,
        method="ANIb",
        program=blastn_tool.exe_path.stem,
        version=blastn_tool.version,
        fragsize=config_anib_args["fragsize"],
        create_db=True,
    )
    # Record the FASTA files in the genomes table _before_ call snakemake
    log_genome(
        database=db,
        fasta=list(check_input_stems(config_anib_args["indir"]).values()),
    )
    assert db.is_file()

    # Run snakemake wrapper
    runner = SnakemakeRunner("snakemake_anib.smk")
    runner.run_workflow(anib_targets_blastn, config_anib_args, workdir=Path(tmp_path))

    # Check output against target fixtures
    for fname in anib_targets_blastn:
        assert fname.suffix == ".tsv", fname
        assert Path(fname).is_file()
        assert (anib_blastn / fname.name).is_file()
        assert filecmp.cmp(fname, anib_blastn / fname.name)

    log_run(
        fasta=config_anib_args["indir"].glob("*.f*"),
        database=db,
        status="Complete",
        name="Test case",
        cmdline="blah blah blah",
        method="ANIb",
        program=blastn_tool.exe_path.stem,
        version=blastn_tool.version,
        fragsize=config_anib_args["fragsize"],
        kmersize=None,
        minmatch=None,
        create_db=False,
    )
    compare_matrices(db, dir_anib_results)
