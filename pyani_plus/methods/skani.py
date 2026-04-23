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


def parse_skani(filename: Path) -> tuple[str, str, float, float, float, str, str]:
    """Parse the skani output files to return a pairwise ANI result.

    :param filename: Path to the skani output file to parse.

    Returns (query path, subject path, ANI percentage, query AF, subject AF, query name, subject name).

    Example skani output file format (tab-delimited):

    ```
    Ref_file	Query_file	ANI	Align_fraction_ref	Align_fraction_query	Ref_name	Query_name
    ../fixtures/viral_example/MGV-GENOME-0264574.fas	../fixtures/viral_example/MGV-GENOME-0264574.fas	\
        100.00	99.65	99.65	MGV_MGV-GENOME-0264574	MGV_MGV-GENOME-0264574
    ```
    """
    # Expected skani header format (tab-delimited):
    expected_header = [
        "Ref_file",
        "Query_file",
        "ANI",
        "Align_fraction_ref",
        "Align_fraction_query",
        "Ref_name",
        "Query_name",
    ]

    # Parse file
    with filename.open() as ifh:
        header = ifh.readline().rstrip("\n").split("\t")
        if header != expected_header:
            msg = f"Unexpected skani output file format in '{filename}'"
            log_sys_exit(logging.getLogger(__name__), msg)
        line = ifh.readline().rstrip("\n").split("\t")
        return (
            line[1],  # query path
            line[0],  # subject path
            float(line[2]),  # ANI percentage
            float(line[4]),  # query AF
            float(line[3]),  # subject AF
            line[6],  # query name
            line[5],  # subject name
        )
