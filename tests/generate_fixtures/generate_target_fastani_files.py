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
"""Generate target files for pyani-plus fastANI tests.

This script can be run with ``./generate_target_fastani_files.py`` in the script's
directory, or from the project root directory via ``make fixtures``. It will
regenerate and potentially modify test input files under the fixtures directory.

This script generates target files for fastani comparisons.
Genomes are compared in both directions (forward and reverse)
using fastANI.
"""

# Imports
import subprocess
from decimal import Decimal
from itertools import product
from pathlib import Path

import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser

from pyani_plus.tools import get_fastani
from pyani_plus.utils import file_md5sum

# Parameters (eg, input sequences, fastANI outputs, k-mer sizes...)
INPUT_DIR = Path("../fixtures/viral_example")
FASTANI_DIR = Path("../fixtures/fastani/targets")
MATRIX_DIR = Path("../fixtures/fastani/matrices")
FRAG_LEN = 3000
KMER_SIZE = 16
MIN_FRAC = 0.2

# Remove pre-existing fixtures before regenerating new ones.
# This is to help with if and when we change the
# example sequences being used.
for file in FASTANI_DIR.glob("*.fastani"):
    file.unlink()

# Running comparisons
inputs = {_.stem: _ for _ in Path(INPUT_DIR).glob("*.f*")}
comparisons = product(inputs, inputs)

fastani = get_fastani()
print(f"Using fastANI {fastani.version} at {fastani.exe_path}")

for genomes in comparisons:
    stem = "_vs_".join(genomes)
    subprocess.run(
        [
            fastani.exe_path,
            "-q",
            inputs[genomes[0]],
            "-r",
            inputs[genomes[1]],
            "-o",
            FASTANI_DIR / (stem + ".fastani"),
            "--fragLen",
            str(FRAG_LEN),
            "-k",
            str(KMER_SIZE),
            "--minFraction",
            str(MIN_FRAC),
        ],
        check=True,
    )

md5dict = {str(file): file_md5sum(file) for file in inputs.values()}
hashes = sorted(set(md5dict.values()))
n = len(hashes)
assert len(inputs) == n, "Duplicate files!"


def fasta_seq_len(fasta_filename: Path) -> int:
    """Get total number of bases/letters in a FASTA file."""
    length = 0
    with Path(fasta_filename).open() as handle:
        for _title, seq in SimpleFastaParser(handle):
            length += len(seq)
    return length


lengths = {file_md5sum(file): fasta_seq_len(file) for file in inputs.values()}


def write_matrix(filename: Path, hashes: list[str], values: np.array) -> None:
    """Write an NxN TSV matrix of values labelled by MD5 hashes."""
    n = len(hashes)
    assert values.shape == (n, n)
    with filename.open("w") as handle:
        handle.write("\t" + "\t".join(hashes) + "\n")
        for row, md5 in enumerate(hashes):
            handle.write(
                md5 + "\t" + "\t".join(str(values[row, col]) for col in range(n)) + "\n"
            )


matrix_ani_string = np.full((n, n), "?", np.dtypes.StringDType)
matrix_orthologous_matches = np.full((n, n), -1, int)
matrix_fragments = np.full((n, n), -1, int)

print("Now parsing the fastANI output to generate expected matrices")
for file in FASTANI_DIR.glob("*.fastani"):
    with file.open() as handle:
        fields = handle.readline().rstrip("\n").split("\t")
        assert len(fields) == 5, f"Bad input {file}"  # noqa: PLR2004
        row = hashes.index(md5dict[fields[0]])  # query
        col = hashes.index(md5dict[fields[1]])  # subject
        # This is to avoid 99.8332 becoming 0.9983329999999999 and so on:
        matrix_ani_string[row, col] = str(Decimal(fields[2]) / 100)
        matrix_orthologous_matches[row, col] = float(fields[3])
        matrix_fragments[row, col] = float(fields[4])

assert matrix_orthologous_matches.min() >= 0
assert matrix_fragments.min() >= 0
print("Now calculating derived values and writing matrices")

# This is an approximation to ANIm style coverage calculated using bp:
matrix_coverage = matrix_orthologous_matches / matrix_fragments

# This is a proxy value, unmatched matrix_fragments rather than bp:
matrix_sim_errors = matrix_fragments - matrix_orthologous_matches

# Element-wise multiplication, must cast string ANI to float (range 0 to 1)
matrix_hadamard = matrix_coverage * matrix_ani_string.astype(float)

# Can't give a very good estimate here:
matrix_aln_lengths = np.full((n, n), 0, int)
for row, query in enumerate(hashes):
    for col, subject in enumerate(hashes):
        matrix_aln_lengths[row, col] = int(
            matrix_coverage[row, col] * min(lengths[query], lengths[subject])
        )

write_matrix(MATRIX_DIR / "matrix_identity.tsv", hashes, matrix_ani_string)
write_matrix(MATRIX_DIR / "matrix_coverage.tsv", hashes, matrix_coverage)
write_matrix(MATRIX_DIR / "matrix_aln_lengths.tsv", hashes, matrix_aln_lengths)
write_matrix(MATRIX_DIR / "matrix_hadamard.tsv", hashes, matrix_hadamard)
write_matrix(MATRIX_DIR / "matrix_sim_errors.tsv", hashes, matrix_sim_errors)

print("Done")
