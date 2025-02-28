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
"""Test snakemake workflow for sourmash.

These tests are intended to be run from the repository root using:

pytest -v or make test
"""

import json
from pathlib import Path

# Required to support pytest automated testing
from pyani_plus.private_cli import log_run, prepare_genomes
from pyani_plus.public_cli import cli_sourmash, resume
from pyani_plus.tools import get_sourmash

from . import compare_db_matrices

KMERSIZE = 31
SCALED = 300  # default scaled=1000 not suitable for the 3 viruses


def compare_sourmash_sig_files(file1: Path, file2: Path) -> bool:
    """Compare two .sig files, ignoring the path part of the filename entry."""
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
                assert Path(entry1[key]).name == Path(entry2[key]).name
            else:
                assert entry1[key] == entry2[key], (
                    f"{key} {entry1[key]!r}!={entry2[key]!r}"
                )

    return True


def test_sketch_prepare(
    input_genomes_tiny: Path,
    tmp_path: str,
) -> None:
    """Test sourmash sketch via the prepare-genomes command."""
    tmp_dir = Path(tmp_path)
    cache = tmp_dir / "cache"
    cache.mkdir()

    tool = get_sourmash()

    tmp_db = tmp_dir / "sig-prepare.db"
    log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus sourmash ...",
        status="Testing",
        name="Testing sourmash prepare-genomes",
        method="sourmash",
        program=tool.exe_path.name,
        version=tool.version,
        kmersize=KMERSIZE,
        extra=f"scaled={SCALED}",
        create_db=True,
    )

    # Run prepare-genomes command...
    prepare_genomes(
        database=tmp_db,
        run_id=1,
        cache=cache,
    )

    # Check output against target fixtures
    for expected in (input_genomes_tiny / "intermediates/sourmash").glob("*.sig"):
        generated = cache / f"sourmash_k={KMERSIZE}_scaled={SCALED}" / expected.name
        assert compare_sourmash_sig_files(expected, generated)

    resume(database=tmp_db, cache=cache)
    compare_db_matrices(tmp_db, input_genomes_tiny / "matrices")


def test_compare_rule_bad_align(
    input_genomes_bad_alignments: Path,
    tmp_path: str,
) -> None:
    """Test sourmash compare snakemake wrapper (bad_alignments).

    Checks that the compare rule in the sourmash snakemake wrapper gives the
    expected output.
    """
    tmp_dir = Path(tmp_path)
    cache = tmp_dir / "cache"
    cache.mkdir()

    # Setup the cache as if prepare-genomes had made the signatures
    cache = cache / f"sourmash_k={KMERSIZE}_scaled={SCALED}"
    cache.mkdir()
    for sig in (input_genomes_bad_alignments / "intermediates/sourmash").glob("*.sig"):
        (cache / sig.name).symlink_to(sig)

    tmp_db = tmp_dir / "bad-align.db"

    cli_sourmash(
        database=tmp_db,
        create_db=True,
        fasta=input_genomes_bad_alignments,
        cache=cache,
        kmersize=KMERSIZE,
        scaled=SCALED,
        temp=tmp_dir,
    )
    compare_db_matrices(tmp_db, input_genomes_bad_alignments / "matrices")
