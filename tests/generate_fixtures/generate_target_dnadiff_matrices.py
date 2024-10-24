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
"""Generate target matrices for pyani-plus dnadiff tests.

This script can be run with ``./generate_target_anim_files.py`` in the
script's directory, or from the project root directory via ``make fixtures``.
It will regenerate and potentially modify test input files under the
fixtures directory.

This script generates target matrices for dnadiff method comparisons from
.report files. First few lines of the .report file look like this:

```
MGV-GENOME-0264574.fas MGV-GENOME-0266457.fna
NUCMER

                               [REF]                [QRY]
[Sequences]
TotalSeqs                          1                    1
AlignedSeqs               1(100.00%)           1(100.00%)
UnalignedSeqs               0(0.00%)             0(0.00%)

[Bases]
TotalBases                     39253                39594
AlignedBases           39169(99.79%)        39176(98.94%)
UnalignedBases             84(0.21%)           418(1.06%)

[Alignments]
1-to-1                             2                    2
TotalLength                    59174                59187
AvgLength                   29587.00             29593.50
AvgIdentity                    99.63                99.63

M-to-M                             2                    2
TotalLength                    59174                59187
AvgLength                   29587.00             29593.50
AvgIdentity                    99.63                99.63
```
Here, we focus on extracting AlignedBases and genome coverage
from line[11], and AvgIdentity from line[23].
"""

import re
from decimal import Decimal
from pathlib import Path

import pandas as pd

from pyani_plus import utils


def parse_dnadiff_report(dnadiff_report: Path) -> tuple[int, Decimal, Decimal]:
    """Return dnadiff values for AlignedBases, genome coverage (QRY) and average identity.

    These coverage and identity are converted to be in the range 0 to 1.
    """
    with Path.open(dnadiff_report) as file:
        lines = file.readlines()

    lines_to_retain = [11, 23]
    lines_of_interest = [lines[i].strip() for i in lines_to_retain]

    aligned_bases = int(
        re.findall(r"(\d+)\s*\(\d+\.\d+%\)\s*$", lines_of_interest[0])[0]
    )
    query_coverage = (
        Decimal(re.findall(r"(\d+\.\d+)%\s*\)$", lines_of_interest[0])[0]) / 100
    )
    avg_identity = Decimal(re.findall(r"(\d+\.\d+)\s*$", lines_of_interest[1])[0]) / 100

    return (aligned_bases, query_coverage, avg_identity)


# Constructing a matrix where the MD5 hashes of test genomes are used as both column names and index.
genome_hashes = {
    file.stem: utils.file_md5sum(file)
    for file in Path("../fixtures/viral_example/").glob("*.f*")
}
sorted_hashes = sorted(genome_hashes.values())
aln_lengths_matrix = pd.DataFrame(index=sorted_hashes, columns=sorted_hashes)
coverage_matrix = pd.DataFrame(index=sorted_hashes, columns=sorted_hashes)
identity_matrix = pd.DataFrame(index=sorted_hashes, columns=sorted_hashes)

# Appending information to matrices
report_files = Path("../fixtures/dnadiff/targets/dnadiff_reports/").glob("*.report")

for file in report_files:
    query, subject = file.stem.split("_vs_")
    query_hash = genome_hashes[query]
    subject_hash = genome_hashes[subject]
    aligned_bases, query_coverage, avg_identity = parse_dnadiff_report(file)
    aln_lengths_matrix.loc[query_hash, subject_hash] = aligned_bases
    coverage_matrix.loc[query_hash, subject_hash] = query_coverage
    identity_matrix.loc[query_hash, subject_hash] = avg_identity

matrices_directory = "../fixtures/dnadiff/matrices/"
Path(matrices_directory).mkdir(parents=True, exist_ok=True)

aln_lengths_matrix.to_csv(matrices_directory + "matrix_aln_lengths.tsv", sep="\t")
coverage_matrix.to_csv(matrices_directory + "matrix_coverage.tsv", sep="\t")
identity_matrix.to_csv(matrices_directory + "matrix_identity.tsv", sep="\t")
