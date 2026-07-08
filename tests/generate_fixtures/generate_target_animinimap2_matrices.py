#!/usr/bin/env python3
#
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
"""Generate target matrices for pyani-plus sourmash tests.

This script can be run with
``./generate_target_sourmash_matrices.py <path_to_inputs_dir> <path_to_output_dir>``
in the script's directory, or from the project root directory via ``make fixtures``.
It will regenerate and potentially modify test input files under the
fixtures directory.

This script generates target matrices for sourmash method comparisons from
sourmash compare .csv files.
"""

import gzip
import subprocess
import sys
import tempfile
from itertools import product
from pathlib import Path

import pandas as pd

from pyani_plus.methods.animinimap2 import parse_minimap2_paf_file
from pyani_plus.tools import get_minimap2
from pyani_plus.utils import fasta_bytes_iterator, file_md5sum, filename_stem

# Paths to directories (eg, input, output)
INPUT_DIR, OUT_DIR = Path(sys.argv[1]), Path(sys.argv[2])

# Constructing a matrix where the stems of test genomes are used as both column names and index.
sorted_stems = sorted(filename_stem(str(file)) for file in INPUT_DIR.glob("*.f*"))
stem_hashes = {
    filename_stem(str(file)): file_md5sum(file) for file in INPUT_DIR.glob("*.f*")
}
aln_lengths_matrix = pd.DataFrame(index=sorted_stems, columns=sorted_stems)
coverage_matrix = pd.DataFrame(index=sorted_stems, columns=sorted_stems)
identity_matrix = pd.DataFrame(index=sorted_stems, columns=sorted_stems)

minimap2 = get_minimap2()
print(f"Using minimap2 version {minimap2.version} at {minimap2.exe_path}")

# Running comparisons (all vs all)
inputs = {filename_stem(str(_)): _ for _ in sorted(Path(INPUT_DIR).glob("*.f*"))}
comparisons = product(inputs, inputs)

# Generate and parse .report files
for query, subject in comparisons:
    in_fasta = inputs[query]
    with (
        gzip.open(in_fasta)
        if in_fasta.suffix == ".gz"
        else in_fasta.open("rb") as handle
    ):
        query_length = sum(len(s) for (_, s) in fasta_bytes_iterator(handle))

    stem = f"{stem_hashes[query]}_vs_{stem_hashes[subject]}"
    # Running minimap2 and saving output to temp file
    with tempfile.TemporaryDirectory() as tmp:
        paf = tmp + "/" + stem + ".paf"
        subprocess.run(
            [
                minimap2.exe_path,
                "-x",
                "asm20",
                "-o",
                paf,
                inputs[subject],
                inputs[query],
            ],
            check=True,
        )
        # Extracting values from PAF TSV file
        query_aligned_bases, subject_aligned_bases, identity = parse_minimap2_paf_file(
            Path(paf)
        )

    # Assign values to matrices
    aln_lengths_matrix.loc[query, subject] = (
        None if query_aligned_bases == 0 else query_aligned_bases
    )
    coverage_matrix.loc[query, subject] = (
        None if query_aligned_bases == 0 else float(query_aligned_bases) / query_length
    )
    identity_matrix.loc[query, subject] = None if identity == 0 else identity

# Save matrices
matrices_directory = OUT_DIR
matrices_directory.mkdir(parents=True, exist_ok=True)

aln_lengths_matrix.to_csv(matrices_directory / "ANIminimap2_aln_lengths.tsv", sep="\t")
coverage_matrix.to_csv(matrices_directory / "ANIminimap2_coverage.tsv", sep="\t")
identity_matrix.to_csv(matrices_directory / "ANIminimap2_identity.tsv", sep="\t")
