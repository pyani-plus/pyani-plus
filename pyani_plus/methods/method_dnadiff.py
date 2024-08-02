# (c) University of Strathclyde 2024-present
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2024-present University of Strathclyde
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
"""Code to implement the dnadiff Average Nucleotide Identity (ANI) method.

This script calculates ANI using the method implemented by the dnadiff
software, available at https://github.com/garviz/MUMmer/blob/master/dnadiff.

All input FASTA format files are compared against each other in both directions
using NUCmer with the `--maxmatch` parameter. To replicate the AlignedBases
value reported by dnadiff, delta-filter (`-m` parameters), show-diff (`-rH`),
and show-coords (default settings) are also run (see snakemake_dnadiff.smk).

The generated files are parsed to obtain values for AlignedBases (alignment
length), average nucleotide identity (ANI) percentages, and aligned percentage
(of the whole genome) for each pairwise comparison.
"""

# Set Up
from pathlib import Path

import pandas as pd  # type: ignore  # noqa: PGH003
from Bio import SeqIO  # type: ignore  # noqa: PGH003


def parse_dnadiff(
    genome_seq: Path,
    mcoords_file: Path,
    rdiff_file: Path,
) -> tuple[int, float, float]:  # type: ignore  # noqa: PGH003
    """Return (reference alignment length, average identity, genome coverage).

    :param genome_seq: Path to inout FASTA files
    :param mcoords_file: Path to the .mcoords file
    :param rdiff_file: Path to the .rdiff file
    """
    # Assigning values
    aligned_bases = 0

    # Step 1. Getting basic sequence information (eg. total number of bases, lengths of individual sequences)
    total_bases = 0

    records = list(SeqIO.parse(genome_seq, "fasta"))
    refs = {record.id: len(record.seq) for record in records}

    for sequence_length in refs.values():
        total_bases += sequence_length

    # Step 2. Check which sequences in the given genome were aligned at least once and update the number of
    # aligned bases by incrementing it with the total length of each unique aligned sequence.
    # Retriving information from M-to-M alignments (mcoords)
    column_names = [
        "S1",
        "E1",
        "S2",
        "E2",
        "LEN 1",
        "LEN 2",
        "% IDY",
        "ref_tag",
        "query_tag",
    ]
    mcoords = pd.read_csv(
        mcoords_file,
        sep=r"\s+\|\s+|\s+",  # Custom separator to handle spaces and pipes
        skiprows=4,  # Skip the first 4 rows
        names=column_names,
        comment="=",  # Ignore the line starting with '=',
        engine="python",
    )

    seen_ref_seq = []
    for _index, row in mcoords.iterrows():
        if row["ref_tag"] not in seen_ref_seq:
            aligned_bases += refs[row["ref_tag"]]
            seen_ref_seq.append(row["ref_tag"])

    # Step 3: Identify lengths of gaps in the aligned sequences and update the AlignedBases values
    # Information retrived from rdiff files
    rdiff = pd.read_csv(
        Path(rdiff_file),
        sep="\t",
        names=[f"col_{_}" for _ in range(7)],
    )

    for _index, row in rdiff.iterrows():
        gap = row["col_4"]
        if row["col_1"] != "DUP" and gap > 0:
            aligned_bases -= gap

    # Step 4: Calculate coverage percentage (aligned bases/total bases * 100%)
    coverage = aligned_bases / total_bases * 100

    # Step 5: Calculate average identity from .mcoords file
    sum_identity = 0
    sum_aligment_lengths = 0
    for _index, row in mcoords.iterrows():
        sum_identity += (row["% IDY"] / 100) * (row["LEN 1"] + row["LEN 2"])
        sum_aligment_lengths += row["LEN 1"] + row["LEN 2"]
    average_identity = round(sum_identity / sum_aligment_lengths * 100, 2)

    return (aligned_bases, coverage, average_identity)
