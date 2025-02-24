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
"""Test snakemake workflow for fastANI.

These tests are intended to be run from the repository root using:

make test
"""

import filecmp
import os
import re
import shutil
from pathlib import Path

import pytest

from pyani_plus.private_cli import log_run
from pyani_plus.tools import get_fastani
from pyani_plus.workflows import (
    ToolExecutor,
    run_snakemake_with_progress_bar,
)

from . import compare_db_matrices


@pytest.fixture
def config_fastani_args(
    fastani_targets_outdir: Path,
    input_genomes_tiny: Path,
    snakemake_cores: int,
    tmp_path: str,
) -> dict:
    """Return configuration settings for testing snakemake fastANI rule."""
    return {
        "db": Path(tmp_path) / "db.slqite",
        "run_id": 1,  # by construction
        "fastani": get_fastani().exe_path,
        "outdir": fastani_targets_outdir,
        "indir": input_genomes_tiny,
        "cores": snakemake_cores,
        "fragsize": 3000,
        "kmersize": 16,
        "minmatch": 0.2,
    }


def compare_fastani_files(file1: Path, file2: Path) -> bool:
    """Compare two fastANI files.

    This function expects two text files as input and returns True if the content
    of the files is the same, and False if the two files differ.

    As the Path to both Query and Reference might be different,
    we will only consider file name.
    """
    with file1.open() as if1, file2.open() as if2:
        # return False if any errors found
        try:
            for line1, line2 in zip(if1.readlines(), if2.readlines(), strict=True):
                if line1.split("\t")[2:] != line2.split("\t")[2:]:
                    return False
        except ValueError:
            return False

    # by default return True
    return True


def test_rule_fastani(
    capsys: pytest.CaptureFixture[str],
    input_genomes_tiny: Path,
    fastani_targets_outdir: Path,
    config_fastani_args: dict,
    tmp_path: str,
) -> None:
    """Test fastANI snakemake wrapper.

    Checks that the fastANI rule in the fastANI snakemake wrapper gives the
    expected output.
    """
    tmp_dir = Path(tmp_path)

    # Assuming this will match but worker nodes might have a different version
    fastani_tool = get_fastani()

    # Setup minimal test DB
    db = config_fastani_args["db"]
    assert not db.is_file()
    log_run(
        fasta=config_fastani_args["indir"],  # i.e. input_genomes_tiny
        database=db,
        status="Testing",
        name="Test case",
        cmdline="blah blah blah",
        method="fastANI",
        program=fastani_tool.exe_path.stem,
        version=fastani_tool.version,
        fragsize=config_fastani_args["fragsize"],
        kmersize=config_fastani_args["kmersize"],
        minmatch=config_fastani_args["minmatch"],
        create_db=True,
    )
    assert db.is_file()
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    # Run snakemake wrapper
    run_snakemake_with_progress_bar(
        executor=ToolExecutor.local,
        workflow_name="compute_column.smk",
        targets=[
            fastani_targets_outdir / f"all_vs_{s.stem}.fastani"
            for s in input_genomes_tiny.glob("*.f*")
        ],
        params=config_fastani_args,
        working_directory=Path(tmp_path),
        temp=Path(tmp_path),
    )

    # Check the intermediate files

    for file in (input_genomes_tiny / "intermediates/fastANI").glob("*_vs_*.fastani"):
        assert filecmp.cmp(file, tmp_dir / file), f"Wrong fastANI output in {file.name}"

    compare_db_matrices(db, input_genomes_tiny / "matrices")


def test_duplicate_stems(
    fastani_targets_outdir: Path,
    config_fastani_args: dict,
    tmp_path: str,
) -> None:
    """Test fastANI snakemake wrapper with duplicated stems in inputs.

    Checks this error condition is caught by duplicating the main fastani
    test, but using a modified copy of the input directory under temp.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(fastani_targets_outdir, ignore_errors=True)

    dup_input_dir = Path(tmp_path) / "duplicated_stems"
    dup_input_dir.mkdir()
    stems = set()
    for sequence in config_fastani_args["indir"].glob("*.f*"):
        # For every input FASTA file, make two versions - XXX.fasta and XXX.fna
        os.symlink(sequence, dup_input_dir / (sequence.stem + ".fasta"))
        os.symlink(sequence, dup_input_dir / (sequence.stem + ".fna"))
        stems.add(sequence.stem)
    dup_config = config_fastani_args.copy()
    dup_config["indir"] = dup_input_dir
    msg = f"Duplicated stems found for {sorted(stems)}. Please investigate."

    generated_targets = [
        fastani_targets_outdir / f"{q}_vs_{s}.fastani"
        for q in dup_input_dir.glob("*.f*")
        for s in dup_input_dir.glob("*.f*")
    ]

    with pytest.raises(ValueError, match=re.escape(msg)):
        run_snakemake_with_progress_bar(
            executor=ToolExecutor.local,
            workflow_name="compute_column.smk",
            targets=generated_targets,
            params=dup_config,
            working_directory=Path(tmp_path),
            temp=Path(tmp_path),
        )
