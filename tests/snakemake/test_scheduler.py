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
"""Tests for the pyani_plus/snakemake/snakemake_scheduler.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

import os
import re
from pathlib import Path

import pytest

from pyani_plus.snakemake import snakemake_scheduler


def test_input_stems(
    input_genomes_small: Path,
) -> None:
    """Verify expected stem:filename mapping from FASTA directory."""
    mapping = snakemake_scheduler.check_input_stems(str(input_genomes_small))
    assert mapping == {
        "NC_002696": input_genomes_small / "NC_002696.fasta",
        "NC_010338": input_genomes_small / "NC_010338.fna",
        "NC_011916": input_genomes_small / "NC_011916.fas",
        "NC_014100": input_genomes_small / "NC_014100.fna",
    }


def test_duplicate_stems(
    input_genomes_small: Path,
    tmp_path: str,
) -> None:
    """Verify expected error message with duplicated FASTA stems."""
    dup_input_dir = Path(tmp_path) / "duplicated_stems"
    dup_input_dir.mkdir()
    stems = set()
    for sequence in input_genomes_small.glob("*.f*"):
        # For every input FASTA file, make two versions - XXX.fasta and XXX.fna
        os.symlink(sequence, dup_input_dir / (sequence.stem + ".fasta"))
        os.symlink(sequence, dup_input_dir / (sequence.stem + ".fna"))
        stems.add(sequence.stem)

    msg = f"Duplicated stems found for {sorted(stems)}. Please investigate."
    with pytest.raises(ValueError, match=re.escape(msg)):
        snakemake_scheduler.check_input_stems(str(dup_input_dir))
