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


class EnumModeSourmash(str, Enum):
    """Enum for the --mode command line argument passed to sourmash."""

    max_containment = "max-containment"  # default
    containment = "containment"


MODE = EnumModeSourmash.max_containment  # constant for CLI default
SCALED = 1000
KMER_SIZE = 31  # default


def parse_sourmash_compare_csv(
    compare_file: Path, filename_to_hash: dict[str, str]
) -> Iterator[tuple[str, str, float | None]]:
    """Parse sourmash compare all-vs-all CSV output.

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
        for row, subject in enumerate(hashes):
            values = handle.readline().rstrip("\n").split(",")
            if values[row] != "1.0":
                # This could happen due to a bug in sourmash, or a glitch in the
                # parser if we're not looking at the matrix element we think we are?
                msg = f"Expected sourmash {subject} vs self to be one, not {values[row]!r}"
                raise ValueError(msg)
            for col, ani in enumerate(values):
                yield hashes[col], subject, None if float(ani) == 0.0 else float(ani)
