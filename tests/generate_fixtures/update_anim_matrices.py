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
"""Update ANIm matrices.

This script updates ANIm matrices by incorporating comparisons that return no
alignments, which would otherwise cause issues with legacy pyANI.
This is executed within ./generate_target_anim_matrices.py and may regenerate or
modify ANIm matrices produced by the legacy pyANI pipeline.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from pyani_plus import utils
from pyani_plus.methods import method_anim

# Load legacy pyANI matrices and add missing rows and columns to accommodate comparisons with no alignments
genome_hashes = {
    file.stem: utils.file_md5sum(file)
    for file in Path("../fixtures/viral_example/").glob("*.f*")
}
sorted_hashes = sorted(genome_hashes.values())
aln_lengths_matrix = pd.read_csv(
    "../fixtures/anim/matrices/matrix_aln_lengths.tsv", sep="\t", index_col=0
).reindex(index=sorted_hashes, columns=sorted_hashes, fill_value=0)
coverage_matrix = pd.read_csv(
    "../fixtures/anim/matrices/matrix_coverage.tsv", sep="\t", index_col=0
).reindex(index=sorted_hashes, columns=sorted_hashes)
identity_matrix = pd.read_csv(
    "../fixtures/anim/matrices/matrix_identity.tsv", sep="\t", index_col=0
).reindex(index=sorted_hashes, columns=sorted_hashes)
sim_errors = pd.read_csv(
    "../fixtures/anim/matrices/matrix_sim_errors.tsv", sep="\t", index_col=0
).reindex(index=sorted_hashes, columns=sorted_hashes, fill_value=0)
hadamard_matrix = pd.read_csv(
    "../fixtures/anim/matrices/matrix_hadamard.tsv", sep="\t", index_col=0
).reindex(index=sorted_hashes, columns=sorted_hashes)

filter_files = Path("../fixtures/anim/targets/filter/").glob("*.filter")

for file in filter_files:
    if "MGV-GENOME-0357962" in file.stem:
        query_length = 87285
        query, subject = file.stem.split("_vs_")
        query_hash = genome_hashes[query]
        subject_hash = genome_hashes[subject]
        query_aligned_bases, subject_aligned_bases, identity, similarity_errors = (
            method_anim.parse_delta(file)
        )
        aln_lengths_matrix.loc[query_hash, subject_hash] = query_aligned_bases
        coverage_matrix.loc[query_hash, subject_hash] = (
            float(query_aligned_bases) / query_length
        )
        identity_matrix.loc[query_hash, subject_hash] = (
            identity if identity is not None else np.nan
        )
        sim_errors.loc[query_hash, subject_hash] = similarity_errors
        hadamard_matrix.loc[query_hash, subject_hash] = (
            identity * (float(query_aligned_bases) / query_length)
            if identity is not None
            else np.nan
        )


matrices_directory = "../fixtures/anim/matrices/"
aln_lengths_matrix.to_csv(matrices_directory + "matrix_aln_lengths.tsv", sep="\t")
coverage_matrix.to_csv(matrices_directory + "matrix_coverage.tsv", sep="\t")
identity_matrix.to_csv(matrices_directory + "matrix_identity.tsv", sep="\t")
sim_errors.to_csv(matrices_directory + "matrix_sim_errors.tsv", sep="\t")
hadamard_matrix.to_csv(matrices_directory + "matrix_hadamard.tsv", sep="\t")
