#!/usr/bin/env python3
#
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
"""Generate target matrices for pyani-plus sourmash tests.

This script can be run with
``./generate_target_sourmash_matrices.py <path_to_inputs_dir> <path_to_output_dir>``
in the script's directory, or from the project root directory via ``make fixtures``.
It will regenerate and potentially modify test input files under the
fixtures directory.

This script generates target matrices for sourmash method comparisons from
sourmash compare .csv files.
"""

import sys
from pathlib import Path

import pandas as pd

from pyani_plus.utils import file_md5sum, filename_stem

# Paths to directories (eg, input, output)
INPUT_DIR, OUT_DIR = Path(sys.argv[1]), Path(sys.argv[2])


def parse_compare_files(compare_file: Path) -> float:
    """Return the top-right value from the sourmash compare .csv file."""
    compare_results = pd.read_csv(
        Path(compare_file),
        sep=",",
        index_col=False,
    )

    return compare_results.iloc[0, -1]


hash_to_stem = {file_md5sum(_): filename_stem(str(_)) for _ in INPUT_DIR.glob("*.f*")}

# Generate target identity matrix for pyani-plus sourmash tests
matrix = pd.read_csv(INPUT_DIR / "intermediates/sourmash/sourmash.csv")
# To match default pyANI-plus output, must map MD5 used in sourmash.cvs
# to filename stems
matrix.columns = [hash_to_stem[_] for _ in matrix]
matrix.index = matrix.columns
matrix = matrix.sort_index(axis=0).sort_index(axis=1).T
# Note we transpose sourmash.csv compared to our own matrices!
#
# Quoting https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-compare-compare-many-signatures
#
#   Note: compare by default produces a symmetric similarity matrix that can be
#   used for clustering in downstream tasks. With --containment, however, this
#   matrix is no longer symmetric and cannot formally be used for clustering.
#
#   The containment matrix is organized such that the value in row A for column B
#   is the containment of the B'th sketch in the A'th sketch, i.e.
#
#   C(A, B) = B.contained_by(A)
#
# See also the sourmash branchwater output which has explicit query and subject
# columns, and query_containment and subject_containment.
#
# I interpret this as column B is the query, row A is the subject - which is the
# transpose of the convention we are using.

# If ANI values can't be estimated (eg. sourmash 0.0) report None instead
matrix = matrix.replace(0.0, None)

matrices_directory = OUT_DIR
Path(matrices_directory).mkdir(parents=True, exist_ok=True)

# We map sourmash query-containment to our query-coverage:
matrix.to_csv(matrices_directory / "sourmash_coverage.tsv", sep="\t")

# Now convert this to max-containment which maps to our identity:
matrix = matrix.where(matrix > matrix.T, matrix.T)
matrix.to_csv(matrices_directory / "sourmash_identity.tsv", sep="\t")
