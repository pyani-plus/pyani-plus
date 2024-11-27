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

from pyani_plus.private_cli import log_run
from pyani_plus.tools import get_blastn, get_makeblastdb
from pyani_plus.workflows import (
    ToolExecutor,
    run_snakemake_with_progress_bar,
)

from . import compare_db_matrices


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
        "run_id": 1,  # by construction
        "blastn": get_blastn().exe_path,
        "makeblastdb": get_makeblastdb().exe_path,
        "outdir": anib_targets_outdir,
        "indir": input_genomes_tiny,
        "cores": snakemake_cores,
        "fragsize": "1020",
    }


def test_snakemake_rule_fragments(
    input_genomes_tiny: Path,
    anib_targets_outdir: Path,
    config_anib_args: dict,
    tmp_path: str,
) -> None:
    """Test overall ANIb snakemake wrapper."""
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(anib_targets_outdir, ignore_errors=True)

    anib_targets_fragments = [
        anib_targets_outdir / (_.stem + "-fragments.fna")
        for _ in input_genomes_tiny.glob("*.f*")
    ]

    # Run snakemake wrapper
    run_snakemake_with_progress_bar(
        executor=ToolExecutor.local,
        workflow_name="snakemake_anib.smk",
        targets=anib_targets_fragments,
        params=config_anib_args,
        working_directory=Path(tmp_path),
    )

    # Check output against target fixtures
    for fname in anib_targets_fragments:
        assert fname.name.endswith("-fragments.fna"), fname
        assert Path(fname).is_file()
        assert filecmp.cmp(
            fname, input_genomes_tiny / "intermediates/ANIb" / fname.name
        )


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
    input_genomes_tiny: Path,
    anib_targets_outdir: Path,
    config_anib_args: dict,
    tmp_path: str,
) -> None:
    """Test overall ANIb snakemake wrapper."""
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(anib_targets_outdir, ignore_errors=True)

    anib_targets_blastdb = [
        anib_targets_outdir / (_.stem + ".njs") for _ in input_genomes_tiny.glob("*.f*")
    ]

    # Run snakemake wrapper
    run_snakemake_with_progress_bar(
        executor=ToolExecutor.local,
        workflow_name="snakemake_anib.smk",
        targets=anib_targets_blastdb,
        params=config_anib_args,
        working_directory=Path(tmp_path),
    )

    # Check output against target fixtures
    for fname in anib_targets_blastdb:
        assert fname.suffix == ".njs", fname
        assert Path(fname).is_file()
        compare_blast_json(
            fname, input_genomes_tiny / "intermediates/ANIb" / fname.name
        )


def test_snakemake_rule_blastn(
    input_genomes_tiny: Path,
    anib_targets_outdir: Path,
    config_anib_args: dict,
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
    log_run(
        fasta=config_anib_args["indir"],  # i.e. input_genomes_tiny
        database=db,
        status="Testing",
        name="Test case",
        cmdline="blah blah blah",
        method="ANIb",
        program=blastn_tool.exe_path.stem,
        version=blastn_tool.version,
        fragsize=config_anib_args["fragsize"],
        create_db=True,
    )
    assert db.is_file()

    expected_targets = list((input_genomes_tiny / "intermediates/ANIb").glob("*.tsv"))
    generated_targets = [anib_targets_outdir / f.name for f in expected_targets]

    # Run snakemake wrapper
    run_snakemake_with_progress_bar(
        executor=ToolExecutor.local,
        workflow_name="snakemake_anib.smk",
        targets=generated_targets,
        params=config_anib_args,
        working_directory=Path(tmp_path),
    )

    # Check output against target fixtures
    for expected, generated in zip(expected_targets, generated_targets, strict=False):
        assert generated.is_file()
        assert filecmp.cmp(expected, generated)
    compare_db_matrices(db, input_genomes_tiny / "matrices")
