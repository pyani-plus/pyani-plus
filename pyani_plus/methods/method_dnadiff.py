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
value reported by dnadiff, delta-filter (`-m` parameters), show-diff (`-qH`),
and show-coords (default settings) are also run (see snakemake_dnadiff.smk).

The generated files are parsed to obtain values for AlignedBases (alignment
length), average nucleotide identity (ANI) percentages, and aligned percentage
(of the whole genome) for each pairwise comparison.
"""

# Set Up
from pathlib import Path
from typing import NamedTuple

import pandas as pd
from Bio import SeqIO

from pyani_plus.methods.method_anim import get_genome_length


class ComparisonResultDnadiff(NamedTuple):
    """A single comparison result for dnadiff."""

    qname: str  # query sequence name
    rname: str  # reference sequence name
    q_aligned_bases_with_gaps: int  # aligned base count of reference sequence with gaps
    avg_id: float  # average nucleotide identity (as a percentage)
    q_length: int  # total base count of query sequence
    r_length: int  # total base count of reference sequence
    alignment_gaps: int
    aligned_bases: int  # tota base count in an alignment
    r_cov: float


def get_record_lengths(genome_seq: Path) -> dict:
    """Return the lengths of records from the given FASTA files,
    keyed by the sequence ID.

    :param genome_seq: Path to the input FASTA file.
    """  # noqa: D205
    records = list(SeqIO.parse(genome_seq, "fasta"))
    return {record.id: len(record.seq) for record in records}


def parse_mcoords(mcoords_file: Path, record_lengths: dict) -> list[float, int]:
    """Parse mcoords file and return avg ID% and number of aligned bases with gaps.

    :parama mcoords_file: Path to the mcoords_file
    :param record_lengths: Path to the input FASTA file.
    """
    mcoords = pd.read_csv(
        Path(mcoords_file),
        sep="\t",
        names=[f"col_{_}" for _ in range(13)],
    )

    sum_identity = 0
    sum_aligment_lengths = 0
    avg_identity = float(0)
    seen_ref_seq = []
    aligned_bases_with_gaps = 0

    # Calculate average nucleotide identity
    for _index, row in mcoords.iterrows():
        sum_identity += (row["col_6"] / 100) * (row["col_4"] + row["col_5"])
        sum_aligment_lengths += row["col_4"] + row["col_5"]
        avg_identity = round(sum_identity / sum_aligment_lengths * 100, 2)

    # Get number of aligned bases with gaps
        if row["col_12"] not in seen_ref_seq:
            aligned_bases_with_gaps += record_lengths[row["col_12"]]
            seen_ref_seq.append(row["col_12"])

    return [avg_identity, aligned_bases_with_gaps]


def parse_qdiff(qdiff_file: Path) -> int:
    """Parse qdiff and return number of gaps in the alignment.

    :param qdiff_file: Path to the qdiff file
    """
    qdiff = pd.read_csv(
        Path(qdiff_file),
        sep="\t",
        names=[f"col_{_}" for _ in range(7)],
    )

    gap_lengths = 0
    for _index, row in qdiff.iterrows():
        gap = row["col_4"]
        if row["col_1"] != "DUP" and gap > 0:
            gap_lengths += gap

    return gap_lengths


def collect_dnadiff_results_directory(
    mcoords_file: Path,
    qdiff_file: Path,
    indir: Path,
) -> ComparisonResultDnadiff:
    """Return a ComparisonResultDnadiff for a completed dnadiff comparison."""
    files = {fasta.stem: fasta for fasta in indir.iterdir() if fasta.is_file()}

    rname, qname = mcoords_file.stem.split("_vs_")
    r_genome_length = get_genome_length(files[rname])
    q_genome_length = get_genome_length(files[qname])
    q_record_lengths = get_record_lengths(files[qname])
    avg_identity, aligned_bases_with_gaps = parse_mcoords(mcoords_file, q_record_lengths)
    gaps = parse_qdiff(qdiff_file)


    return ComparisonResultDnadiff(
        rname=rname,
        qname=qname,
        avg_id=avg_identity,
        q_aligned_bases_with_gaps=aligned_bases_with_gaps,
        r_length=r_genome_length,
        q_length=q_genome_length,
        alignment_gaps=gaps,
        aligned_bases=aligned_bases_with_gaps - gaps,
        r_cov=(aligned_bases_with_gaps - gaps) / q_genome_length * 100,
    )
