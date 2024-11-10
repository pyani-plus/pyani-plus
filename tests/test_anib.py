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
"""Tests for the ANIb implementation.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, tools
from pyani_plus.methods import method_anib

from . import get_matrix_entry


def test_bad_path(tmp_path: str) -> None:
    """Confirm giving an empty path etc fails."""
    with pytest.raises(
        FileNotFoundError, match="No such file or directory: '/does/not/exist'"
    ):
        method_anib.fragment_fasta_files([Path("/does/not/exist")], Path(tmp_path))
    with pytest.raises(
        FileNotFoundError, match="No such file or directory: '/does/not/exist'"
    ):
        method_anib.parse_blastn_file(Path("/does/not/exist"))


def test_empty_path(tmp_path: str) -> None:
    """Confirm giving an empty path etc fails."""
    with pytest.raises(ValueError, match="No sequences found in /dev/null"):
        method_anib.fragment_fasta_files([Path("/dev/null")], Path(tmp_path))
    # Note it is valid to have an empty BLASTN TSV file (there are no headers)


def test_parse_blastn_empty() -> None:
    """Check parsing of empty BLASTN tabular file."""
    assert method_anib.parse_blastn_file(Path("/dev/null")) == (0, 0, 0)


def test_parse_blastn_bad(input_genomes_tiny: Path) -> None:
    """Check parsing something which isn't a BLAST TSV files fails."""
    with pytest.raises(
        ValueError, match="Found 1 columns in .*/MGV-GENOME-0264574.fas, not 15"
    ):
        method_anib.parse_blastn_file(input_genomes_tiny / "MGV-GENOME-0264574.fas")


def test_parse_blastn_bad_query(tmp_path: str) -> None:
    """Check parsing TSV not using frag#### fails."""
    fake = Path(tmp_path) / "fake.tsv"
    with fake.open("w") as handle:
        handle.write("\t".join(["bad-query", "subject"] + ["0"] * 13))
    with pytest.raises(
        ValueError,
        match="BLAST output should be using fragmented queries, not bad-query",
    ):
        method_anib.parse_blastn_file(fake)


def test_parse_blastn(anib_blastn: Path) -> None:
    """Check parsing of BLASTN tabular file."""
    # Function returns tuple of mean percentage identity, total alignment length, and
    # total mismatches/gaps:
    assert method_anib.parse_blastn_file(
        anib_blastn / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.tsv"
    ) == (0.9945938461538462, 39169, 215)
    # $ md5sum tests/fixtures/viral_example/*.f*
    # 689d3fd6881db36b5e08329cf23cecdd tests/fixtures/viral_example/MGV-GENOME-0264574.fas
    # 78975d5144a1cd12e98898d573cf6536 tests/fixtures/viral_example/MGV-GENOME-0266457.fna
    # 5584c7029328dc48d33f95f0a78f7e57 tests/fixtures/viral_example/OP073605.fasta
    #
    # Example is query 689d3fd6881db36b5e08329cf23cecdd vs subject 78975d5144a1cd12e98898d573cf6536
    #
    # Expected matrices from pyani v0.2 give us expected values of:
    # identity 0.9945938461538463, aln_length 39169, sim_errors 215
    # Note final digit wobble from 3 (old) to 2 (new)


def test_missing_db(tmp_path: str, input_genomes_tiny: Path, anib_blastn: Path) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.log_anib(
            database=tmp_db,
            run_id=1,
            # These are for the comparison table
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            blastn=anib_blastn / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.tsv",
        )


def test_bad_query_or_subject(
    tmp_path: str, input_genomes_tiny: Path, anib_blastn: Path
) -> None:
    """Mismatch between query or subject FASTA in fastANI output and commandline."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --query-fasta .*/MGV-GENOME-0266457.fna"
            " but query in blastn filename was MGV-GENOME-0264574"
        ),
    ):
        private_cli.log_anib(
            database=tmp_db,
            run_id=1,
            # These are for the comparison table
            query_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            blastn=anib_blastn / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.tsv",
        )

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --subject-fasta .*/MGV-GENOME-0264574.fas"
            " but subject in blastn filename was MGV-GENOME-0266457"
        ),
    ):
        private_cli.log_anib(
            database=tmp_db,
            run_id=1,
            # These are for the comparison table
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            blastn=anib_blastn / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.tsv",
        )


def test_logging_anib(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
    anib_blastn: Path,
    dir_anib_results: Path,
) -> None:
    """Check can log a ANIb comparison to DB."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_blastn()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus anib ...",
        status="Testing",
        name="Testing logging_anib",
        method="ANIb",
        program=tool.exe_path.stem,
        version=tool.version,
        fragsize=method_anib.FRAGSIZE,  # will be used by default in log_anib
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    private_cli.log_anib(
        database=tmp_db,
        run_id=1,
        # These are for the comparison table
        query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
        subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
        blastn=anib_blastn / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.tsv",
    )

    # Check the recorded comparison values
    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 1
    comp = session.query(db_orm.Comparison).one()
    query = "689d3fd6881db36b5e08329cf23cecdd"  # MGV-GENOME-0264574.fas
    subject = "78975d5144a1cd12e98898d573cf6536"  # MGV-GENOME-0266457.fna
    pytest.approx(
        comp.identity,
        get_matrix_entry(dir_anib_results / "matrix_identity.tsv", query, subject),
    )
    pytest.approx(
        comp.aln_length,
        get_matrix_entry(dir_anib_results / "matrix_aln_lengths.tsv", query, subject),
    )
    pytest.approx(
        comp.sim_errors,
        get_matrix_entry(dir_anib_results / "matrix_sim_errors.tsv", query, subject),
    )
    pytest.approx(
        comp.cov_query,
        get_matrix_entry(dir_anib_results / "matrix_coverage.tsv", query, subject),
    )
    session.close()
    tmp_db.unlink()
