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

from pathlib import Path


def parse_fastani_file(filename: Path) -> tuple[Path, Path, float, int, int]:
    """Parse a single-line fastANI output file extracting key fields as a tuple.

    Return (ref genome, query genome, ANI estimate, orthologous matches,
    sequence fragments) tuple.

    :param filename: Path, path to the input file

    Extracts the ANI estimate (which we return in the range 0 to 1), the
    number of orthologous matches (int), and the number of sequence
    fragments considered from the fastANI output file (int).

    We assume that all fastANI comparisons are pairwise: one query and
    one reference file. The fastANI file should contain a single line.

    fastANI *can* produce multi-line output, if a list of query/reference
    files is given to it.
    """
    with filename.open() as handle:
        line = handle.readline().strip().split()
    if not line:  # No file content; either run failed or no detectable similarity
        msg = f"Input file {filename} is empty"
        raise ValueError(msg)
    return (
        Path(line[0]),
        Path(line[1]),
        0.01 * float(line[2]),
        int(line[3]),
        int(line[4]),
    )
