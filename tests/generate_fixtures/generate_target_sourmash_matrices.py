#!/usr/bin/env python3
#
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

from pyani_plus import utils

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


# Constructing a matrix where the MD5 hashes of test genomes are used as both column names and index.
genome_hashes = {file.stem: utils.file_md5sum(file) for file in INPUT_DIR.glob("*.f*")}
sorted_hashes = sorted(genome_hashes.values())

# Generate target identity matrix for pyani-plus sourmash tests
identity_matrix = pd.read_csv(INPUT_DIR / "intermediates/sourmash/sourmash.csv")
identity_matrix.columns = [genome_hashes[Path(_).stem] for _ in identity_matrix]
identity_matrix.index = identity_matrix.columns

# If ANI values can't be estimated (eg. sourmash 0.0) report None instead
identity_matrix = identity_matrix.replace(0.0, None)


matrices_directory = OUT_DIR
Path(matrices_directory).mkdir(parents=True, exist_ok=True)
identity_matrix.to_csv(matrices_directory / "sourmash_identity.tsv", sep="\t")
