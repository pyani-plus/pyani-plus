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
"""Tests for the fastANI implementation.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, tools
from pyani_plus.methods.method_fastani import parse_fastani_file

from . import get_matrix_entry


def test_bad_path() -> None:
    """Confirm giving an empty path etc fails."""
    with pytest.raises(FileNotFoundError, match="No such file or directory:"):
        parse_fastani_file(Path("/does/not/exist"))


def test_empty_path() -> None:
    """Confirm giving an empty path etc fails."""
    with pytest.raises(ValueError, match="Input file /dev/null is empty"):
        parse_fastani_file(Path("/dev/null"))


def test_bad_query_or_subject(
    input_genomes_tiny: Path, fastani_targets_indir: Path
) -> None:
    """Mismatch between query or subject FASTA in fastANI output and commandline."""
    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --query-fasta .*/MGV-GENOME-0266457.fna"
            " but query in fastANI file was .*/MGV-GENOME-0264574.fas"
        ),
    ):
        private_cli.log_fastani(
            database=":memory:",
            # These are for the comparison table
            query_fasta=input_genomes_tiny
            / "MGV-GENOME-0266457.fna",  # should be MGV-GENOME-0264574.fas
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            fastani=fastani_targets_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.fastani",
            # These are all for the configuration table:
            fragsize=1000,
            kmersize=51,
            minmatch=0.9,
        )

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --subject-fasta .*/MGV-GENOME-0264574.fas"
            " but query in fastANI file was .*/MGV-GENOME-0266457.fna"
        ),
    ):
        private_cli.log_fastani(
            database=":memory:",
            # These are for the comparison table
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny
            / "MGV-GENOME-0264574.fas",  # should be MGV-GENOME-0266457.fna",
            fastani=fastani_targets_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.fastani",
            # These are all for the configuration table:
            fragsize=1000,
            kmersize=51,
            minmatch=0.9,
        )


def test_missing_db(
    tmp_path: str, input_genomes_tiny: Path, fastani_targets_indir: Path
) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.log_fastani(
            database=tmp_db,
            # These are for the comparison table
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            fastani=fastani_targets_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.fastani",
            # These are all for the configuration table:
            fragsize=1000,
            kmersize=51,
            minmatch=0.9,
        )


def test_logging_fastani(
    tmp_path: str,
    input_genomes_tiny: Path,
    fastani_targets_indir: Path,
    fastani_matrices: Path,
) -> None:
    """Check can log a fastANI comparison to DB."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_fastani()

    private_cli.log_configuration(
        database=tmp_db,
        method="fastANI",
        program=tool.exe_path.stem,
        version=tool.version,
        fragsize=1000,
        kmersize=51,
        minmatch=0.9,
        create_db=True,
    )
    private_cli.log_genome(
        database=tmp_db,
        fasta=[
            input_genomes_tiny / "MGV-GENOME-0264574.fas",
            input_genomes_tiny / "MGV-GENOME-0266457.fna",
        ],
    )

    private_cli.log_fastani(
        database=tmp_db,
        # These are for the comparison table
        query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
        subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
        fastani=fastani_targets_indir
        / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.fastani",
        # These are all for the configuration table:
        fragsize=1000,
        kmersize=51,
        minmatch=0.9,
    )

    # Check the recorded comparison values
    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 1
    comp = session.query(db_orm.Comparison).one()
    query = "689d3fd6881db36b5e08329cf23cecdd"  # MGV-GENOME-0264574.fas
    subject = "78975d5144a1cd12e98898d573cf6536"  # MGV-GENOME-0266457.fna
    # Using approx to avoid 0.9950140000000001 != 0.995014
    pytest.approx(
        comp.identity,
        get_matrix_entry(fastani_matrices / "matrix_identity.tsv", query, subject),
    )
    pytest.approx(
        comp.aln_length,
        get_matrix_entry(fastani_matrices / "matrix_aln_lengths.tsv", query, subject),
    )
    pytest.approx(
        comp.sim_errors,
        get_matrix_entry(fastani_matrices / "matrix_sim_errors.tsv", query, subject),
    )
    pytest.approx(
        comp.cov_query,
        get_matrix_entry(fastani_matrices / "matrix_coverage.tsv", query, subject),
    )
    session.close()
    tmp_db.unlink()
