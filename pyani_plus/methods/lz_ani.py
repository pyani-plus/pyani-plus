# The MIT License
#
# Copyright (c) 2024-2026 University of Strathclyde
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
"""Code to implement the skani Average Nucleotide Identity (ANI) method."""

import logging
from pathlib import Path

from pyani_plus import log_sys_exit


def parse_lzani(
    filename: Path,
) -> tuple[
    tuple[str, str, float | None, float | None],
    tuple[str, str, float | None, float | None],
]:
    """Parse the lz-ani output files to return a pairwise ANI result.

    :param filename: Path to the lz-ani output file to parse.

    As lz-ani carries out forward and reverse comparisons, we return a tuple of two tuples,
    one for the forward comparison and one for the reverse comparison. The internal tuples are:

    (query file, subject file, ANI percentage, query coverage)

    Example lz-ani output file format (tab-delimited):

    ```
    qidx	ridx	query	reference	tani	gani	ani	qcov	num_alns	len_ratio
    1	0	NC_011916.fas.gz	NC_014100.fna.gz	0.512457	0.547429	0.868086	0.630616	7849	0.8684
    0	1	NC_014100.fna.gz	NC_011916.fas.gz	0.512457	0.482088	0.866453	0.556392	8415	0.8684
    ```
    """
    # Expected lz-ani header format (tab-delimited):
    expected_header = [
        "qidx",
        "ridx",
        "query",
        "reference",
        "tani",
        "gani",
        "ani",
        "qcov",
        "num_alns",
        "len_ratio",
    ]

    # Parse file
    with filename.open() as ifh:
        header = ifh.readline().rstrip("\n").split("\t")
        if header != expected_header:
            msg = f"Unexpected lz-ani output file format in '{filename}'"
            log_sys_exit(logging.getLogger(__name__), msg)
        fline = ifh.readline().rstrip("\n").split("\t")  # forward comparison line
        rline = ifh.readline().rstrip("\n").split("\t")  # reverse comparison line

        # Return the pairwise alignment contents
        return (
            ("", "", None, None)
            if fline == [""]
            else (fline[2], fline[3], float(fline[6]), float(fline[7])),
            ("", "", None, None)
            if rline == [""]
            else (rline[2], rline[3], float(rline[6]), float(rline[7])),
        )
