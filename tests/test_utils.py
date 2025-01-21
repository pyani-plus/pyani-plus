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
"""Tests for the pyani_plus/utils.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pathlib import Path

import pytest

from pyani_plus import utils


def test_md5_str() -> None:
    """Confirm our MD5 function works with a gzip filename as a string."""
    assert (
        utils.file_md5sum("tests/fixtures/bacterial_example/NC_002696.fasta.gz")
        == "f19cb07198a41a4406a22b2f57a6b5e7"
    )


def test_md5_path() -> None:
    """Confirm our MD5 function works with a filename as a Path."""
    assert (
        utils.file_md5sum(Path("tests/fixtures/bacterial_example/NC_002696.fasta.gz"))
        == "f19cb07198a41a4406a22b2f57a6b5e7"
    )


def test_md5_invalid() -> None:
    """Confirm our MD5 function failure mode with non-existent filename."""
    with pytest.raises(
        ValueError, match="Input file /does/not/exist.txt is not a file or symlink"
    ):
        utils.file_md5sum("/does/not/exist.txt")


def test_check_output() -> None:
    """Confirm our subprocess wrapper catches expected failures."""
    # I wanted to check the full stderr, but couldn't get it to work.
    # After and outside the context manager capsys.readouterr().out was empty.
    # Inside the context manager, the code never ran...
    with pytest.raises(
        SystemExit,
        match=r'ERROR: Return code 1 from: blastn -task blast\nError: Argument "task". Illegal value',
    ):
        utils.check_output(["blastn", "-task", "blast"])


def test_stage_file(tmp_path: Path) -> None:
    """Check error conditions when staging a file."""
    tmp_dir = Path(tmp_path)

    with pytest.raises(
        SystemExit,
        match=r"ERROR: Missing input file /does/not/exist.in",
    ):
        utils.stage_file(Path("/does/not/exist.in"), Path("/does/not/exist.out"))

    tmp_inp = tmp_dir / "example.fasta.gz"
    tmp_inp.touch()
    tmp_out = tmp_dir / "example.fasta"
    tmp_out.touch()

    with pytest.raises(
        SystemExit,
        match=r"ERROR: Intermediate file .*/example.fasta already exists!",
    ):
        utils.stage_file(tmp_inp, tmp_out)
