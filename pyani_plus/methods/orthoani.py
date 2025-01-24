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
"""Code to implement the orthoANI average nucleotide identity method.

Calculates ANI by the orthoANI method, as described in Lee et al. (2016)
Int J Syst Evol Micr 66: 1100-1103. https://doi.org/10.1099/ijsem.0.000760.

The orthoANI algorithm calculates ANI in three main steps. First, both the
query and reference genomes are fragmented into segments of 1020 bp.
Fragments shorter than 1020 bp are excluded. Next, NCBI BLAST is used
to compare the fragmented query sequences against a database of fragmented
reference sequences. Finally, the ANI value is determined based on reciprocal
best hits between the query and reference genomes.
"""

import gzip
from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser

FRAGSIZE = 1020  # Default orthiANI fragment size


def fragment_fasta_file(
    filename: Path,
    fragmented_fasta: Path,
    method: str,
    fragsize: int = FRAGSIZE,
) -> None:
    """Fragment FASTA file into subsequences of up to the given size.

    For ANIb, any remainder is taken as is (1 ≤ length < fragsize),
    while for orthoANI, the remainders are discarded.

    Accepts gzipped files as input.
    """
    with (
        (
            gzip.open(filename, "rt") if filename.suffix == ".gz" else filename.open()
        ) as in_handle,
        fragmented_fasta.open("w") as out_handle,
    ):
        count = 0
        for title, seq in SimpleFastaParser(in_handle):
            index = 0
            while index < len(seq):
                fragment = seq[index : index + fragsize]
                # Exclude fragments shorter than 1020bp if the method is orthoANI
                if method == "orthoANI" and len(fragment) < fragsize:
                    index += fragsize
                    continue
                count += 1
                out_handle.write(f">frag{count:05d} {title}\n")
                # Line wrap at 60 chars
                for i in range(0, len(fragment), 60):
                    out_handle.write(fragment[i : i + 60] + "\n")
                index += fragsize

    if not count:
        msg = f"No sequences found in {filename}"
        raise ValueError(msg)
