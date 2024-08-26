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
