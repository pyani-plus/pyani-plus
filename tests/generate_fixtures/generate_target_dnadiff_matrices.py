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
"""Generate target matrices for pyani-plus dnadiff tests.

This script can be run with
``./generate_target_anim_files.py <path_to_input_dir> <path_to_output_dir>``
in the script's directory, or from the project root directory via ``make fixtures``.
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
import subprocess
import sys
import tempfile
from decimal import Decimal
from itertools import product
from pathlib import Path

import pandas as pd

from pyani_plus.tools import get_dnadiff
from pyani_plus.utils import file_md5sum

INPUT_DIR, OUT_DIR = Path(sys.argv[1]), Path(sys.argv[2])

# Constructing a matrix where the stems of test genomes are used as both column names and index.
sorted_stems = sorted(file.stem for file in INPUT_DIR.glob("*.f*"))
stem_hashes = {file.stem: file_md5sum(file) for file in INPUT_DIR.glob("*.f*")}
aln_lengths_matrix = pd.DataFrame(index=sorted_stems, columns=sorted_stems)
coverage_matrix = pd.DataFrame(index=sorted_stems, columns=sorted_stems)
identity_matrix = pd.DataFrame(index=sorted_stems, columns=sorted_stems)
sim_errors = pd.DataFrame(index=sorted_stems, columns=sorted_stems)


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


dnadiff = get_dnadiff()
print(f"Using nucmer {dnadiff.version} at {dnadiff.exe_path}")


# Running comparisons (all vs all)
inputs = {_.stem: _ for _ in sorted(Path(INPUT_DIR).glob("*.f*"))}
comparisons = product(inputs, inputs)

# Generate and parse .report files
for query, subject in comparisons:
    stem = f"{stem_hashes[query]}_vs_{stem_hashes[subject]}"
    # Running dnadiff and saving output to temp file
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(
            [
                dnadiff.exe_path,
                "-p",
                tmp + "/" + stem,
                inputs[subject],
                inputs[query],
            ],
            check=True,
        )
        # Extracting values from .report files
        aligned_bases, query_coverage, avg_identity = parse_dnadiff_report(
            Path(tmp + "/" + stem + ".report")
        )

        # Set all values to None if all metrics are zero
        values = (
            (None, None, None, None)
            if aligned_bases == 0 and query_coverage == 0 and avg_identity == 0
            else (
                aligned_bases,
                query_coverage,
                avg_identity,
                round(aligned_bases * (1 - avg_identity)),
            )
        )

        # Assign values to matrices
        (
            aln_lengths_matrix.loc[query, subject],
            coverage_matrix.loc[query, subject],
            identity_matrix.loc[query, subject],
            sim_errors.loc[query, subject],
        ) = values


# Save matrices
matrices_directory = OUT_DIR
matrices_directory.mkdir(parents=True, exist_ok=True)

aln_lengths_matrix.to_csv(matrices_directory / "dnadiff_aln_lengths.tsv", sep="\t")
coverage_matrix.to_csv(matrices_directory / "dnadiff_coverage.tsv", sep="\t")
identity_matrix.to_csv(matrices_directory / "dnadiff_identity.tsv", sep="\t")
sim_errors.to_csv(matrices_directory / "dnadiff_sim_errors.tsv", sep="\t")
