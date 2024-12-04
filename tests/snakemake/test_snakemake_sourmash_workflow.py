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

import json
import shutil  # We need this for filesystem operations
from pathlib import Path

import pandas as pd

# Required to support pytest automated testing
import pytest

from pyani_plus.private_cli import log_run
from pyani_plus.tools import get_sourmash
from pyani_plus.workflows import (
    ToolExecutor,
    run_snakemake_with_progress_bar,
)

from . import compare_db_matrices


@pytest.fixture
def config_sourmash_args(
    snakemake_cores: int,
    tmp_path: str,
) -> dict:
    """Return configuration settings for testing snakemake sourmash rules."""
    return {
        "db": Path(tmp_path) / "db.sqlite",
        "run_id": 1,  # by construction
        # "outdir": ... is dynamic
        # "indir": ... is dynamic
        "cores": snakemake_cores,
        "kmersize": 31,
        "mode": "containment",  # default is max-containment
        "extra": "scaled=300",  # default scaled=1000 not suitable for the 3 viruses
    }


def compare_sourmash_sig_files(file1: Path, file2: Path) -> bool:
    """Compare two .sig files, considering only the stem in the filename field.

    Also ignores name field if missing in either file (branchwater quirk).
    """
    with Path.open(file1) as f1:
        data1 = json.load(f1)

    with Path.open(file2) as f2:
        data2 = json.load(f2)

    assert isinstance(data1, list)
    assert isinstance(data2, list)
    assert len(data1) == len(data2)

    for entry1, entry2 in zip(data1, data2, strict=False):
        assert isinstance(entry1, dict)
        assert isinstance(entry2, dict)
        keys = set(entry1).union(entry2)
        for key in keys:
            if key == "filename":
                assert Path(entry1[key]).stem == Path(entry2[key]).stem
            elif key == "name" and (key not in entry2 or key not in entry1):
                # Known to be missing in branchwater
                pass
            else:
                assert (
                    entry1[key] == entry2[key]
                ), f"{key} {entry1[key]!r}!={entry2[key]!r}"

    return True


def compare_sourmash_ani_files(data1: Path, data2: Path) -> bool:
    """Compare two .csv files returned by sourmash compare."""
    # Read the .csv files into DataFrames
    csv1 = pd.read_csv(data1)
    csv2 = pd.read_csv(data2)

    # Use regex to replace everything before the last '/' with an empty string
    csv1.columns = csv1.columns.str.replace(r".*/", "", regex=True)
    csv2.columns = csv2.columns.str.replace(r".*/", "", regex=True)

    return csv1.equals(csv2)


def compare_sourmash_ani_files_bad_alignments(data1: Path, data2: Path) -> bool:
    """Compare two .csv files returned by sourmash compare."""
    # Read the .csv files into DataFrames
    csv1 = pd.read_csv(data1)
    csv2 = pd.read_csv(data2)

    # Use regex to replace everything before the last '/' with an empty string
    csv1.columns = csv1.columns.str.replace(r".*/", "", regex=True)
    csv2.columns = csv2.columns.str.replace(r".*/", "", regex=True)

    # Note: In sourmash compare, the absence of ANI value estimation is
    # represented as 0.0 for distantly related genomes. In pyANI-plus,
    # these 0.0 values are converted to None to represent NULL values.
    # Ensuring this difference is accounted for when comparing files.
    csv2 = csv2.replace(0.0, None).astype(object)
    csv1 = csv1.replace(0.0, None).astype(object)

    return csv1.equals(csv2)


def test_sketch_rule(
    input_genomes_tiny: Path,
    sourmash_targets_signature_outdir: Path,
    config_sourmash_args: dict,
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
    config["indir"] = input_genomes_tiny

    expected_sigs = list((input_genomes_tiny / "intermediates/sourmash").glob("*.sig"))
    targets = [
        sourmash_targets_signature_outdir / fname.name for fname in expected_sigs
    ]

    # Run snakemake wrapper
    run_snakemake_with_progress_bar(
        executor=ToolExecutor.local,
        workflow_name="snakemake_sourmash.smk",
        targets=targets,
        params=config,
        working_directory=Path(tmp_path),
    )

    # Check output against target fixtures
    for expected, generated in zip(expected_sigs, targets, strict=False):
        assert compare_sourmash_sig_files(expected, generated)


def test_compare_rule(
    capsys: pytest.CaptureFixture[str],
    sourmash_targets_compare_outdir: Path,
    config_sourmash_args: dict,
    input_genomes_tiny: Path,
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
    config["indir"] = input_genomes_tiny

    # Assuming this will match but worker nodes might have a different version
    sourmash_tool = get_sourmash()

    # Setup minimal test DB
    db = config["db"]
    assert not db.is_file()
    log_run(
        fasta=config["indir"],  # i.e. input_genomes_tiny
        database=db,
        cmdline="pyani-plus sourmash ...",
        status="Testing",
        name="Testing sourmash",
        method="sourmash",
        program=sourmash_tool.exe_path.stem,
        version=sourmash_tool.version,
        mode=config["mode"],
        kmersize=config["kmersize"],
        extra=config["extra"],
        create_db=True,
    )
    assert db.is_file()
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    # Run snakemake wrapper
    run_snakemake_with_progress_bar(
        executor=ToolExecutor.local,
        workflow_name="snakemake_sourmash.smk",
        targets=[sourmash_targets_compare_outdir / "sourmash.csv"],
        params=config,
        working_directory=Path(tmp_path),
    )

    # Check output against target fixture
    assert compare_sourmash_ani_files(
        sourmash_targets_compare_outdir / "sourmash.csv",
        input_genomes_tiny / "intermediates/sourmash/sourmash.csv",
    )

    compare_db_matrices(db, input_genomes_tiny / "matrices")


def test_compare_rule_bad_alignments(
    capsys: pytest.CaptureFixture[str],
    sourmash_targets_compare_outdir: Path,
    config_sourmash_args: dict,
    input_genomes_bad_alignments: Path,
    tmp_path: str,
) -> None:
    """Test sourmash compare snakemake wrapper (bad_alignments).

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
    config["indir"] = input_genomes_bad_alignments

    # Assuming this will match but worker nodes might have a different version
    sourmash_tool = get_sourmash()

    # Setup minimal test DB
    db = config["db"]
    assert not db.is_file()
    log_run(
        fasta=config["indir"],  # i.e. input_genomes_tiny
        database=db,
        cmdline="pyani-plus sourmash ...",
        status="Testing",
        name="Testing sourmash",
        method="sourmash",
        program=sourmash_tool.exe_path.stem,
        version=sourmash_tool.version,
        mode=config["mode"],
        kmersize=config["kmersize"],
        extra=config["extra"],
        create_db=True,
    )
    assert db.is_file()
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    # Run snakemake wrapper
    run_snakemake_with_progress_bar(
        executor=ToolExecutor.local,
        workflow_name="snakemake_sourmash.smk",
        targets=[sourmash_targets_compare_outdir / "sourmash.csv"],
        params=config,
        working_directory=Path(tmp_path),
    )

    # Check output against target fixture
    assert compare_sourmash_ani_files_bad_alignments(
        sourmash_targets_compare_outdir / "sourmash.csv",
        input_genomes_bad_alignments / "intermediates/sourmash/sourmash.csv",
    )

    compare_db_matrices(db, input_genomes_bad_alignments / "matrices")


def test_branchwater_sketch_rule(
    input_genomes_tiny: Path,
    config_sourmash_args: dict,
    tmp_path: str,
) -> None:
    """Test sourmash branchwater sketch snakemake wrapper."""
    tmp_dir = Path(tmp_path)

    config = config_sourmash_args.copy()
    config["outdir"] = tmp_dir / "output"
    config["indir"] = input_genomes_tiny

    expected_sigs = list((input_genomes_tiny / "intermediates/sourmash").glob("*.sig"))
    targets = [tmp_dir / "output" / fname.name for fname in expected_sigs]

    # Run snakemake wrapper
    run_snakemake_with_progress_bar(
        executor=ToolExecutor.local,
        workflow_name="snakemake_branchwater.smk",
        targets=targets,
        params=config,
        working_directory=Path(tmp_path),
    )

    # Check output against target fixtures
    for expected, generated in zip(expected_sigs, targets, strict=False):
        assert compare_sourmash_sig_files(expected, generated)


def test_branchwater_compare_rule(
    capsys: pytest.CaptureFixture[str],
    config_sourmash_args: dict,
    input_genomes_tiny: Path,
    tmp_path: str,
) -> None:
    """Test sourmash branchwater compare snakemake wrapper."""
    tmp_dir = Path(tmp_path)

    config = config_sourmash_args.copy()
    config["outdir"] = tmp_dir / "output"
    config["indir"] = input_genomes_tiny

    # Assuming this will match but worker nodes might have a different version
    sourmash_tool = get_sourmash()

    # Setup minimal test DB
    db = config["db"]
    assert not db.is_file()
    log_run(
        fasta=config["indir"],  # i.e. input_genomes_tiny
        database=db,
        cmdline="pyani-plus branchwater ...",
        status="Testing",
        name="Testing branchwater",
        method="branchwater",
        program=sourmash_tool.exe_path.stem,
        version=sourmash_tool.version,
        # not implemented yet, mode=config["mode"],
        kmersize=config["kmersize"],
        extra=config["extra"],
        mode=config["mode"],
        create_db=True,
    )
    assert db.is_file()
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    # Run snakemake wrapper
    run_snakemake_with_progress_bar(
        executor=ToolExecutor.local,
        workflow_name="snakemake_branchwater.smk",
        targets=[tmp_dir / "output/branchwater.csv"],
        params=config,
        working_directory=Path(tmp_path),
    )

    compare_db_matrices(db, input_genomes_tiny / "matrices")
