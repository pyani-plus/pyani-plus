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
        "scaled": 300,
        "kmer": 31,
    }


def compare_sourmash_sig_files(
    file1: Path, file2: Path, key_to_exclude: str = "filename"
) -> bool:
    """Compare two .sig files, excluding filename ifnormation.

    Return True if they are the same or False if they are different.
    """

    def remove_key(data: dict, key_to_remove: str) -> dict:
        """Remove the specified key (file) from a JSON-like structure (dict or list).

        This is because sourmash provides the full path for input genomes in .sig files,
        which will differ when comparing target files with those generated during workflow run.
        """
        if isinstance(data, dict):  # If the data is a dictionary
            # Remove the key if it exists in the dictionary
            data.pop(key_to_remove, None)
            # Call remove_key on each value in the dictionary
            for value in data.values():
                remove_key(value, key_to_remove)
        elif isinstance(data, list):  # If the data is a list
            # Call remove_key on each item in the list
            for item in data:
                remove_key(item, key_to_remove)

        return data

    # Load sketch data from the first .sig file
    with Path.open(file1) as f1:
        data1 = json.load(f1)

    # Load sketch data from the second .sig file
    with Path.open(file2) as f2:
        data2 = json.load(f2)

    # Remove the specified key from both loaded sig files
    cleaned_data1 = remove_key(data1, key_to_exclude)
    cleaned_data2 = remove_key(data2, key_to_exclude)

    # Compare the cleaned .sig data and return True if they are identical, False if different
    return cleaned_data1 == cleaned_data2


def compare_sourmash_ani_files(data1: Path, data2: Path) -> bool:
    """Compare two .csv files returned by sourmash compare."""
    # Read the .csv files into DataFrames
    csv1 = pd.read_csv(data1)
    csv2 = pd.read_csv(data2)

    # Use regex to replace everything before the last '/' with an empty string
    csv1.columns = csv1.columns.str.replace(r".*/", "", regex=True)
    csv2.columns = csv2.columns.str.replace(r".*/", "", regex=True)

    return csv1.equals(csv2)


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
    for fname in sourmash_targets_sig:
        assert Path(fname).suffix == ".sig", fname
        assert Path(fname).is_file()
        assert compare_sourmash_sig_files(
            sourmash_targets_signature_outdir / Path(fname).name,
            sourmash_targets_sig_indir / Path(fname).name,
        )


def test_snakemake_compare_rule(
    sourmash_targets_compare_outdir: Path,
    sourmash_target_ani: list[str],
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
    runner.run_workflow(sourmash_target_ani, config, workdir=Path(tmp_path))

    # Check output against target fixtures
    for fname in sourmash_target_ani:
        assert Path(fname).suffix == ".csv", fname
        assert Path(fname).is_file()
        assert compare_sourmash_ani_files(
            sourmash_targets_compare_outdir / Path(fname).name,
            sourmash_targets_compare_indir / Path(fname).name,
        )
