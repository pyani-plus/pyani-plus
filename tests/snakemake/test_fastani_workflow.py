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
from pathlib import Path

import pytest

from pyani_plus import setup_logger
from pyani_plus.db_orm import connect_to_db
from pyani_plus.private_cli import import_json_comparisons, log_run
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
    logger = setup_logger(None)
    # Run snakemake wrapper
    json_targets = [
        fastani_targets_outdir / f"fastani.run_1.column_{_ + 1}.json" for _ in range(3)
    ]
    run_snakemake_with_progress_bar(
        logger,
        executor=ToolExecutor.local,
        workflow_name="compute_column.smk",
        targets=json_targets,
        database=db,
        working_directory=Path(tmp_path),
        temp=Path(tmp_path),
    )

    # Check the intermediate files
    for file in (input_genomes_tiny / "intermediates/fastANI").glob("*_vs_*.fastani"):
        assert filecmp.cmp(file, tmp_dir / file), f"Wrong fastANI output in {file.name}"

    session = connect_to_db(logger, db)
    for json in json_targets:
        import_json_comparisons(logger, session, json)
    session.close()

    compare_db_matrices(db, input_genomes_tiny / "matrices")
