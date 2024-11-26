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

import subprocess
import tempfile
from decimal import Decimal
from pathlib import Path

import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser

from pyani_plus.tools import get_fastani
from pyani_plus.utils import file_md5sum

# Parameters (eg, input sequences, fastANI outputs, k-mer sizes...)
INPUT_DIR = Path("../fixtures/viral_example")
FASTANI_DIR = Path("../fixtures/fastani/targets")
MATRIX_DIR = Path("../fixtures/viral_example/matrices")
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

fastani = get_fastani()
print(f"Using fastANI {fastani.version} at {fastani.exe_path}")

# Generate values for each row of the matrix
with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as genome_list:
    for genome_filename in inputs.values():
        genome_list.write(f"{genome_filename}\n")
    genome_list.close()
    for stem, genome_filename in inputs.items():
        subprocess.run(
            [
                fastani.exe_path,
                "--ql",
                genome_list.name,
                "-r",
                genome_filename,
                "-o",
                FASTANI_DIR / f"all_vs_{stem}.fastani",
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
        for line in handle:
            fields = line.rstrip("\n").split("\t")
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

# Another estimate:
matrix_aln_lengths = FRAG_LEN * matrix_orthologous_matches

write_matrix(MATRIX_DIR / "fastANI_identity.tsv", hashes, matrix_ani_string)
write_matrix(MATRIX_DIR / "fastANI_coverage.tsv", hashes, matrix_coverage)
write_matrix(MATRIX_DIR / "fastANI_aln_lengths.tsv", hashes, matrix_aln_lengths)
write_matrix(MATRIX_DIR / "fastANI_hadamard.tsv", hashes, matrix_hadamard)
write_matrix(MATRIX_DIR / "fastANI_sim_errors.tsv", hashes, matrix_sim_errors)

print("Done")
