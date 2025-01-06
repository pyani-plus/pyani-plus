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
"""pyani-plus test module."""

from pathlib import Path


def get_matrix_entry(tsv_filename: Path, query_hash: str, subject_hash: str) -> float:
    """Parse the given TSV matrix and pull out a single value as a float."""
    msg = f"Could not find [row={query_hash}, col={subject_hash}] in {tsv_filename}"
    with tsv_filename.open() as handle:
        fields = handle.readline().rstrip("\n").split("\t")
        try:
            col_index = fields.index(subject_hash)
        except ValueError:
            raise ValueError(msg) from None
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == query_hash:
                return float(fields[col_index])
    raise ValueError(msg) from None
