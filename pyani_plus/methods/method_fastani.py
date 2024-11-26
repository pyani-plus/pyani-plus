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
"""Code to wrap the fastANI average nucleotide identity method."""

from collections.abc import Iterator
from pathlib import Path

KMER_SIZE = 16  # Default fastANI k-mer size (max 16)
FRAG_LEN = 3000  # Default fastANI fragment length, mapped to our --fragsize
MIN_FRACTION = 0.2  # Default fastANI min fraction
MAX_RATIO_DIFF = 10.0  # Default fastANI maximum ratio difference


def parse_fastani_file(
    filename: Path, filename_to_hash: dict[str, str]
) -> Iterator[tuple[str, str, float, int, int]]:
    """Parse a multi-line fastANI output file extracting key fields as a tuple.

    Returns tuples of (query genome (hash), reference/subject genome (hash),
    ANI estimate (float), orthologous matches (int), sequence fragments).

    :param filename: Path, path to the input file

    Extracts the ANI estimate (which we return in the range 0 to 1), the
    number of orthologous matches (int), and the number of sequence
    fragments considered from the fastANI output file (int).

    Example fastANI comparison, three queries vs one reference:

    >>> mapping = {
    ...     "MGV-GENOME-0264574.fas": "A",
    ...     "MGV-GENOME-0266457.fna": "B",
    ...     "OP073605.fasta": "C",
    ... }
    >>> fname = Path("tests/fixtures/fastani/targets/all_vs_OP073605.fastani")
    >>> for query, subject, ani, matches, frags in parse_fastani_file(fname, mapping):
    ...     print(f"{query} vs {subject} gave {100*ani:0.1f}%")
    A vs C gave 99.9%
    B vs C gave 99.5%
    C vs C gave 100.0%

    Empty files trigger an exception:

    >>> fname = Path("/dev/null")
    >>> for query, subject, ani, matches, frags in parse_fastani_file(fname, mapping):
    ...     print(f"{query} vs {subject} gave {100*ani:0.1f}%")
    Traceback (most recent call last):
      ...
    ValueError: Input file /dev/null is empty
    """
    empty = True
    with filename.open() as handle:
        for line in handle:
            parts = line.strip().split()
            yield (
                filename_to_hash[Path(parts[0]).name],
                filename_to_hash[Path(parts[1]).name],
                0.01 * float(parts[2]),
                int(parts[3]),
                int(parts[4]),
            )
            empty = False
    if empty:
        msg = f"Input file {filename} is empty"
        raise ValueError(msg)
