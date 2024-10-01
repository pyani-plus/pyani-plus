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

from pyani_plus.private_cli import (
    log_fastani,
    parse_fastani_file,  # may move to methods?
)


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
        log_fastani(
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
            create_db=False,
        )

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --subject-fasta .*/MGV-GENOME-0264574.fas"
            " but query in fastANI file was .*/MGV-GENOME-0266457.fna"
        ),
    ):
        log_fastani(
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
            create_db=False,
        )


def test_missing_db(
    tmp_path: str, input_genomes_tiny: Path, fastani_targets_indir: Path
) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist, but not using --create-db"):
        log_fastani(
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
            create_db=False,
        )

    # This should work:
    log_fastani(
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
        create_db=True,
    )
    tmp_db.unlink()
