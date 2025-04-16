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

from pathlib import Path

# Required to support pytest automated testing
import pytest

from pyani_plus import setup_logger
from pyani_plus.db_orm import connect_to_db
from pyani_plus.private_cli import import_json_comparisons, log_run
from pyani_plus.tools import get_delta_filter, get_nucmer
from pyani_plus.workflows import (
    ToolExecutor,
    run_snakemake_with_progress_bar,
)

from . import compare_db_matrices


@pytest.fixture
def config_anim_args(
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
        # "indir": ... is dynamic
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


def test_rule_ANIm(  # noqa: N802
    input_genomes_tiny: Path,
    config_anim_args: dict,
    anim_targets_outdir: Path,
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
    tmp_dir = Path(tmp_path)

    # Assuming this will match but worker nodes might have a different version
    nucmer_tool = get_nucmer()

    config = config_anim_args.copy()
    config["outdir"] = tmp_path
    config["indir"] = input_genomes_tiny

    # Setup minimal test DB
    db = config_anim_args["db"]
    assert not db.is_file()
    log_run(
        fasta=config["indir"],  # i.e. input_genomes_tiny
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
    logger = setup_logger(None)
    # Run snakemake wrapper
    targets = [
        anim_targets_outdir / f"anim.run_1.column_{_ + 1}.json" for _ in range(3)
    ]
    run_snakemake_with_progress_bar(
        logger,
        executor=ToolExecutor.local,
        workflow_name="compute_column.smk",
        database=db,
        run_id=1,
        targets=targets,
        working_directory=tmp_dir,
        temp=tmp_dir,
    )
    session = connect_to_db(logger, db)
    for json in targets:
        import_json_comparisons(logger, session, json)
    session.close()

    # Check the intermediate files

    # Check delta-filter output against target fixtures
    for fname in (input_genomes_tiny / "intermediates/ANIm").glob("*.filter"):
        generated = next(tmp_dir.glob(f"*/{fname.name}"))
        assert compare_files_with_skip(fname, generated)

    # Check nucmer output against target fixtures
    for fname in (input_genomes_tiny / "intermediates/ANIm").glob("*.delta"):
        generated = next(tmp_dir.glob(f"*/{fname.name}"))
        assert compare_files_with_skip(fname, generated)

    compare_db_matrices(db, input_genomes_tiny / "matrices")


def test_rule_ANIm_bad_align(  # noqa: N802
    input_genomes_bad_alignments: Path,
    config_anim_args: dict,
    anim_targets_outdir: Path,
    tmp_path: str,
) -> None:
    """Test rule ANIm (bad alignments).

    Checks that the ANIm rule in the ANIm snakemake wrapper gives the
    expected output.

    If the output directory exists (i.e. the make clean_tests rule has not
    been run), the tests will automatically pass as snakemake will not
    attempt to re-run the rule. That would prevent us from seeing any
    introduced bugs, so we force re-running the rule by deleting the
    output directory before running the tests.
    """
    tmp_dir = Path(tmp_path)

    # Assuming this will match but worker nodes might have a different version
    nucmer_tool = get_nucmer()

    config = config_anim_args.copy()
    config["outdir"] = tmp_path
    config["indir"] = input_genomes_bad_alignments

    # Setup minimal test DB
    db = config_anim_args["db"]
    assert not db.is_file()
    log_run(
        fasta=config["indir"],  # i.e. input_genomes_bad_alignments
        database=db,
        status="Testing",
        name="Test case",
        cmdline="pyani-plus anib --database ... blah blah blah",
        method="ANIm",
        program=nucmer_tool.exe_path.stem,
        version=nucmer_tool.version,
        mode=config["mode"],
        create_db=True,
    )
    assert db.is_file()
    logger = setup_logger(None)
    # Run snakemake wrapper
    targets = [
        anim_targets_outdir / f"anim.run_1.column_{_ + 1}.json" for _ in range(2)
    ]
    run_snakemake_with_progress_bar(
        logger,
        executor=ToolExecutor.local,
        workflow_name="compute_column.smk",
        database=db,
        run_id=1,
        targets=targets,
        working_directory=tmp_dir,
        temp=tmp_dir,
    )
    session = connect_to_db(logger, db)
    for json in targets:
        import_json_comparisons(logger, session, json)
    session.close()

    # Check delta-filter output against target fixtures
    for fname in (input_genomes_bad_alignments / "intermediates/ANIm").glob("*.filter"):
        generated = next(tmp_dir.glob(f"*/{fname.name}"))
        assert compare_files_with_skip(fname, generated)

    # Check nucmer output against target fixtures
    for fname in (input_genomes_bad_alignments / "intermediates/ANIm").glob("*.delta"):
        generated = next(tmp_dir.glob(f"*/{fname.name}"))
        assert compare_files_with_skip(fname, generated)

    compare_db_matrices(db, input_genomes_bad_alignments / "matrices")
