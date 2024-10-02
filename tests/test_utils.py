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
"""Tests for the pyani_plus/utils.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pathlib import Path

import pytest

from pyani_plus import utils


def test_md5_str() -> None:
    """Confirm our MD5 function works with a filename as a string."""
    assert (
        utils.file_md5sum("tests/fixtures/bacterial_example/NC_002696.fasta")
        == "f19cb07198a41a4406a22b2f57a6b5e7"
    )


def test_md5_path() -> None:
    """Confirm our MD5 function works with a filename as a Path."""
    assert (
        utils.file_md5sum(Path("tests/fixtures/bacterial_example/NC_002696.fasta"))
        == "f19cb07198a41a4406a22b2f57a6b5e7"
    )


def test_md5_invalid() -> None:
    """Confirm our MD5 function failure mode with non-existent filename."""
    with pytest.raises(
        ValueError, match="Input file /does/not/exist.txt is not a file or symlink"
    ):
        utils.file_md5sum("/does/not/exist.txt")
