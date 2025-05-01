# The MIT License
#
# Copyright (c) 2024-2025 University of Strathclyde
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
from pathlib import Path

# Required to support pytest automated testing
import pytest

from pyani_plus import setup_logger
from pyani_plus.db_orm import connect_to_db
from pyani_plus.private_cli import import_json_comparisons, log_run
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


def test_rule_anib(
    input_genomes_tiny: Path,
    anib_targets_outdir: Path,
    config_anib_args: dict,
    tmp_path: str,
) -> None:
    """Test blastn (overall) ANIb snakemake wrapper."""
    tmp_dir = Path(tmp_path)

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
        cmdline="pyani-plus anib blah blah blah",
        method="ANIb",
        program=blastn_tool.exe_path.stem,
        version=blastn_tool.version,
        fragsize=config_anib_args["fragsize"],
        create_db=True,
    )
    assert db.is_file()

    # Run snakemake wrapper
    logger = setup_logger(None)
    json_targets = [
        anib_targets_outdir / f"anib.run_1.column_{_ + 1}.json" for _ in range(3)
    ]
    run_snakemake_with_progress_bar(
        logger,
        executor=ToolExecutor.local,
        workflow_name="compute_column.smk",
        targets=json_targets,
        database=db,
        working_directory=tmp_dir,
        temp=tmp_dir,
    )

    # Check the intermediate files

    for file in (input_genomes_tiny / "intermediates/ANIb").glob("*.f*"):
        assert filecmp.cmp(file, tmp_dir / file), f"Wrong fragmented FASTA {file.name}"

    for file in (input_genomes_tiny / "intermediates/ANIb").glob("*.njs"):
        assert compare_blast_json(file, tmp_dir / file), f"Wrong BLAST DB {file.name}"

    for file in (input_genomes_tiny / "intermediates/ANIb").glob("*_vs_*.tsv"):
        assert filecmp.cmp(file, tmp_dir / file), f"Wrong blastn output in {file.name}"

    session = connect_to_db(logger, db)
    for json in json_targets:
        import_json_comparisons(logger, session, json)
    session.close()

    # Check output against target fixtures
    compare_db_matrices(db, input_genomes_tiny / "matrices")
