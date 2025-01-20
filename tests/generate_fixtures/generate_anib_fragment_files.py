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
"""Generate target files for pyani-plus ANIb (blast) tests.

This script can be run with ``./generate_target_anib_files.py`` in the script's
directory, or from the project root directory via ``make fixtures``. It will
regenerate and potentially modify test input files under the fixtures directory.

This script generates fragmented FASTA files which are used as the query files
with NCBI BLAST+ command blastn against database of the unfragmented FASTA files.
"""

from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser

from pyani_plus.utils import file_md5sum

# Paths to directories (eg, input sequences, delta and filter)
INPUT_DIR = Path("../fixtures/viral_example")
FRAG_DIR = Path("../fixtures/viral_example/intermediates/ANIb")
FRAGSIZE = 1020  # Default ANIb fragment size

# Remove pre-existing fixtures before regenerating new ones.
# This is to help with if and when we change the
# example sequences being used.
for file in FRAG_DIR.glob("*.fna"):
    file.unlink()

# Note flexible on input *.fna vs *.fa vs *.fasta, but fixed on output
for fasta in INPUT_DIR.glob("*.f*"):
    md5 = file_md5sum(fasta)
    output = FRAG_DIR / (md5 + "-fragments.fna")
    count = 0
    with output.open("w") as out_handle, fasta.open("r") as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            index = 0
            while index < len(seq):
                count += 1
                fragment = seq[index : index + FRAGSIZE]
                out_handle.write(f">frag{count:05d} {title}\n")
                # Now line wrap at 60 chars
                for i in range(0, len(fragment), 60):
                    out_handle.write(fragment[i : i + 60] + "\n")
                index += FRAGSIZE
    print(f"Wrote {count} fragments from {fasta.stem}")
