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
"""Code to implement the ANIb average nucleotide identity method.

Calculates ANI by the ANIb method, as described in Goris et al. (2007)
Int J Syst Evol Micr 57: 81-91. https://doi.org/10.1099/ijs.0.64483-0.

The method is based on NCBI BLAST comparisons of fragmented query sequences
against databases of (unfragmented) reference sequences. Thus in an all-vs-all
ANIb comparison of N genomes, we must initially prepare N fragmented genomes,
and separately N BLAST databases of the original genomes, and then nucleotide
BLAST for the N^2 pairwise combinations. This is done with here three snakemake
rules.
"""

from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser

FRAGSIZE = 1020  # Default ANIb fragment size


def fragment_fasta_files(
    fasta: list[Path], outdir: Path, fragsize: int = FRAGSIZE
) -> list[Path]:
    """Fragment FASTA files into subsequences of up to the given size.

    The output files are named ``<stem>-fragments.fna`` regardless of the
    input file extension (typically ``.fna``, ``.fa`` or ``.fasta``).

    Any remainder is taken as is (1 <= length < fragsize).

    Returns a list of the output fragmented FASTA files.
    """
    fragmented_files = []
    for filename in fasta:
        if not isinstance(filename, Path):
            msg = f"Expected a Path object in list of FASTA files, got {filename!r}"
            raise TypeError(msg)
        if not filename.is_file():
            msg = f"Cannot fragment {filename!s}, file not found"
            raise ValueError(msg)
        frag_filename = outdir / (filename.stem + "-fragments.fna")
        with filename.open() as in_handle, frag_filename.open("w") as out_handle:
            count = 0
            for title, seq in SimpleFastaParser(in_handle):
                index = 0
                while index < len(seq):
                    count += 1
                    fragment = seq[index : index + fragsize]
                    out_handle.write(f">frag{count:05d} {title}\n")
                    # Now line wrap at 60 chars
                    for i in range(0, len(fragment), 60):
                        out_handle.write(fragment[i : i + 60] + "\n")
                    index += fragsize
        fragmented_files.append(frag_filename)
    return fragmented_files
