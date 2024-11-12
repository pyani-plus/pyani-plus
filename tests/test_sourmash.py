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
"""Tests for the sourmash implementation.

These tests are intended to be run from the repository root using:

make test
"""

# Required to support pytest automated testing
from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, tools
from pyani_plus.methods import method_sourmash

from . import get_matrix_entry


def test_missing_db(tmp_path: str, sourmash_targets_compare_indir: Path) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.log_sourmash(
            database=tmp_db,
            run_id=1,
            compare=sourmash_targets_compare_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv",
        )


def test_bad_query_or_subject(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
    sourmash_targets_compare_indir: Path,
) -> None:
    """Mismatch between query or subject FASTA in sourmash compare output and commandline.

    Defining a run with two files, but then parsing compare output with different FASTA.
    """
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    tmp_input = Path(tmp_path) / "input"
    tmp_input.mkdir()
    for filename in input_genomes_tiny.glob("MGV-*.f*"):
        alias = tmp_input / (filename.name)
        alias.symlink_to(filename)

    tool = tools.get_sourmash()

    private_cli.log_run(
        fasta=tmp_input,
        database=tmp_db,
        cmdline="pyani-plus sourmash ...",
        status="Testing",
        name="Testing log_sourmash",
        method="sourmash",
        program=tool.exe_path.name,
        version=tool.version,
        mode=method_sourmash.MODE,
        kmersize=method_sourmash.KMER_SIZE,
        extra="scaled=" + str(method_sourmash.SCALED),
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    with pytest.raises(
        SystemExit,
        match=(
            "CSV file .*/MGV-GENOME-0264574_vs_OP073605.csv contained reference"
            " to 'OP073605.fasta' which is not in the run"
        ),
    ):
        private_cli.log_sourmash(
            database=tmp_db,
            run_id=1,
            compare=sourmash_targets_compare_indir
            / "MGV-GENOME-0264574_vs_OP073605.csv",
        )


def test_not_all_vs_all(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
    sourmash_targets_compare_indir: Path,
) -> None:
    """Mismatch where CSV is only a subset.

    Defining a run with 3 files, but then parsing compare output only 2 FASTA.
    """
    tmp_db = Path(tmp_path) / "three.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_sourmash()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus sourmash ...",
        status="Testing",
        name="Testing log_sourmash",
        method="sourmash",
        program=tool.exe_path.name,
        version=tool.version,
        mode=method_sourmash.MODE,
        kmersize=method_sourmash.KMER_SIZE,
        extra="scaled=" + str(method_sourmash.SCALED),
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    with pytest.raises(
        SystemExit,
        match=("Expected sourmash compare CSV to have 3 columns, not 2"),
    ):
        private_cli.log_sourmash(
            database=tmp_db,
            run_id=1,
            compare=sourmash_targets_compare_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv",
        )


def test_logging_wrong_version(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
    sourmash_targets_compare_indir: Path,
) -> None:
    """Check mismatched sourmash version fails."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus sourmash ...",
        status="Testing",
        name="Testing log_sourmash",
        method="sourmash",
        program="sourmash",
        version="42",
        mode=method_sourmash.MODE,
        kmersize=method_sourmash.KMER_SIZE,
        extra="scaled=" + str(method_sourmash.SCALED),
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    with pytest.raises(
        SystemExit,
        match="ERROR: Run configuration was sourmash 42 but we have sourmash 4.",
    ):
        private_cli.log_sourmash(
            database=tmp_db,
            run_id=1,
            compare=sourmash_targets_compare_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv",
        )


def test_logging_sourmash(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
    sourmash_targets_compare_indir: Path,
    dir_sourmash_results: Path,
) -> None:
    """Check can log a sourmash comparison to DB."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    # Mock input folder of just 2 files
    tmp_input = Path(tmp_path) / "input"
    tmp_input.mkdir()
    for filename in input_genomes_tiny.glob("MGV-*.f*"):
        alias = tmp_input / (filename.name)
        alias.symlink_to(filename)

    tool = tools.get_sourmash()

    private_cli.log_run(
        fasta=tmp_input,
        database=tmp_db,
        cmdline="pyani-plus sourmash ...",
        status="Testing",
        name="Testing log_sourmash",
        method="sourmash",
        program=tool.exe_path.stem,
        version=tool.version,
        mode=method_sourmash.MODE,
        kmersize=method_sourmash.KMER_SIZE,
        extra="scaled=" + str(method_sourmash.SCALED),
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    private_cli.log_sourmash(
        database=tmp_db,
        run_id=1,
        compare=sourmash_targets_compare_indir
        / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv",
    )

    # Check the recorded comparison values
    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 4  # noqa: PLR2004
    # Should we just compare as a matrix?
    for comp in session.query(db_orm.Comparison):
        pytest.approx(
            comp.identity,
            get_matrix_entry(
                dir_sourmash_results / "matrix_identity.tsv",
                comp.query_hash,
                comp.subject_hash,
            ),
        )
    session.close()
    tmp_db.unlink()
