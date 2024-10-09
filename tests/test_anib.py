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

from pyani_plus.methods import method_anib


def test_bad_path(tmp_path: str) -> None:
    """Confirm giving an empty path etc fails."""
    with pytest.raises(
        FileNotFoundError, match="No such file or directory: '/does/not/exist'"
    ):
        method_anib.fragment_fasta_files([Path("/does/not/exist")], Path(tmp_path))


def test_empty_path(tmp_path: str) -> None:
    """Confirm giving an empty path etc fails."""
    with pytest.raises(ValueError, match="No sequences found in /dev/null"):
        method_anib.fragment_fasta_files([Path("/dev/null")], Path(tmp_path))
    # Note it is valid to have an empty BLASTN TSV file (there are no headers)
