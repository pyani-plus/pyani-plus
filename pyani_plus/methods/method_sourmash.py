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
"""Code to implement the sourmash Average Nucleotide Identity (ANI) method."""

# Set Up
import sys
from collections.abc import Iterator
from enum import Enum
from pathlib import Path

import pandas as pd


class EnumModeSourmash(str, Enum):
    """Enum for the --mode command line argument passed to sourmash."""

    max_containment = "max-containment"  # default
    containment = "containment"


MODE = EnumModeSourmash.max_containment  # constant for CLI default
SCALED = 1000
KMER_SIZE = 31  # default


def parse_sourmash_compare_csv(
    compare_file: Path, filename_to_hash: dict[str, str]
) -> Iterator[tuple[str, str, float]]:
    """Parse sourmash copare all-vs-all CSV output.

    Returns tuples of (query_hash, subject_hash, estimated ANI).
    """
    with compare_file.open() as handle:
        headers = handle.readline().rstrip("\n").split(",")
        if len(headers) != len(filename_to_hash):
            msg = (
                "Expected sourmash compare CSV to have"
                f" {len(filename_to_hash)} columns, not {len(headers)}"
            )
            sys.exit(msg)
        try:
            hashes = [filename_to_hash[Path(_).name] for _ in headers]
        except KeyError as err:
            msg = f"CSV file {compare_file} contained reference to {err!s} which is not in the run"
            sys.exit(msg)
        for row, query in enumerate(hashes):
            values = handle.readline().rstrip("\n").split(",")
            if values[row] != "1.0":
                msg = f"Expected sourmash {query} vs self to be one, not {values[row]}"
                raise ValueError(msg)
            for col, ani in enumerate(values):
                yield query, hashes[col], float(ani)


def parse_compare(compare_file: Path) -> float:
    """Parse sourmash compare .csv output extracting estimated ANI as float.

    :param compare_file: Path to the compare_file.

    Extracts the ANI estimate (which we return in the range 0 to 1).

    Assumes all sourmash comparisons are pairwise, involving one query and
    one reference file. The sourmash compare file should contain a DataFrame
    with two rows and two columns, e.g.:

        genome_1.fna    genome_2.fas
        1.000000        0.997901
        0.997901        1.000000

    To extract the estimated ANI value, we call `df.iloc[0, -1]`
    to obtain the top right corner value.

    Since sourmash *can* produce multi-comparison output when given a list
    of .sig files, it may be necessary to parse the file by extracting
    the upper triangle of values from the sourmash compare output using
    the `np.triu` function.
    """
    compare_results = pd.read_csv(
        Path(compare_file),
        sep=",",
        index_col=False,
    )

    return compare_results.iloc[0, -1]
