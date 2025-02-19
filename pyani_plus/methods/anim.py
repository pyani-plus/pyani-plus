# Based on MIT licensed code from pyANI,
# Copyright (c) 2016-2019 The James Hutton Institute
# Copyright (c) 2019-2024 University of Strathclyde
#
# The MIT License
#
# Copyright (c) 2019-2025 University of Strathclyde
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
#
"""Code to implement the ANIm average nucleotide identity method.

Calculates ANI by the ANIm method, as described in Richter et al (2009)
Proc Natl Acad Sci USA 106: 19126-19131 doi:10.1073/pnas.0906412106.

All input FASTA format files are compared against each other, pairwise,
using NUCmer (binary location must be provided). NUCmer output will be stored
in a specified output directory.

The NUCmer .delta file output is parsed to obtain an alignment length
and similarity error count for every unique region alignment. These are
processed to give matrices of aligned sequence lengths, similarity error
counts, average nucleotide identity (ANI) percentages, and minimum aligned
percentage (of whole genome) for each pairwise comparison.
"""

from collections import defaultdict
from pathlib import Path

import intervaltree  # type: ignore  # noqa: PGH003

from pyani_plus.public_cli_args import EnumModeANIm

MODE = EnumModeANIm.mum  # constant for CLI default


def get_aligned_bases_count(aligned_regions: dict) -> int:
    """Return count of aligned bases across a set of aligned regions.

    Returns the count of bases in the input sequence that participate in
    an alignment. This is calculated using an intervaltree to merge overlapping
    regions and sum the total number of bases in the merged regions.

    :param aligned_regions: dict of aligned regions
    """
    aligned_bases = 0
    for region in aligned_regions.values():
        tree = intervaltree.IntervalTree.from_tuples(region)
        tree.merge_overlaps(strict=False)
        for interval in tree:
            aligned_bases += interval.end - interval.begin + 1

    return aligned_bases


def parse_delta(
    filename: Path,
) -> tuple[int, int, float | None, int] | tuple[None, None, None, None]:
    """Return (reference alignment length, query alignment length, average identity, similarity errors).

    :param filename: Path to the input .delta file

    Calculates similarity errors and the aligned lengths for reference
    and query and average nucleotide identity, and returns the cumulative
    total for each as a tuple.

    The delta file format contains seven numbers in the lines of interest:
    see http://mummer.sourceforge.net/manual/ for specification

    - start on query
    - end on query
    - start on target
    - end on target
    - error count (non-identical, plus indels)
    - similarity errors (non-positive match scores)
        [NOTE: with PROmer this is equal to error count]
    - stop codons (always zero for nucmer)

    We report ANIm identity by finding an average across all alignments using
    the following formula:

    sum of weighted identical bases / sum of aligned bases from each fragment

    For example:

    reference.fasta query.fasta
    NUCMER
    >ref_seq_A ref_seq_B 40 40
    1 10 1 11 5 5 0
    -1
    0
    15 20 25 30 0 0 0

    The delta file tells us there are two alignments. The first alignment runs from base 1
    to base 10 in the reference sequence, and from base 1 to 11 in the query sequence
    with a similarity error of 5. The second alignment runs from base 15 to 20 in
    the reference, and base 25 to 30 in the query with 0 similarity errors. To calculate
    the %ID, we can:

    - Find the number of all aligned bases from each sequence:
    aligned reference bases region 1 = 10 - 1 + 1 = 10
    aligned query bases region 1 = 11 - 1 + 1 = 11
    aligned reference bases region 2 = 20 - 15 + 1 = 6
    aligned query bases region 2 = 30 - 25 + 1 = 6

    - Find weighted identical bases
    alignment 1 identity weighted = (10 + 11) - (2 * 5) = 11
    alignment 2 identity weighted = (6 + 6) - (2 * 0) = 12

    - Calculate %ID
    (11 + 12) / (10 + 11 + 6 + 6) = 0.696969696969697

    To calculate alignment lengths, we extract the regions of each alignment
    (either for query or reference) provided in the .delta file and merge the overlapping
    regions with IntervalTree. Then, we calculate the total sum of all aligned regions.
    """
    sim_error = 0  # Hold a count for similarity errors

    regions_ref = defaultdict(list)  # Hold a dictionary for query regions
    regions_qry = defaultdict(list)  # Hold a dictionary for query regions

    aligned_bases = 0  # Hold a count of aligned bases for each sequence
    weighted_identical_bases = 0  # Hold a count of weighted identical bases

    # Flag to check if any ">" line exists
    # This is to capture comparisons that result in no alignment
    contains_alignment = False

    # Ideally we wouldn't read the whole file into memory at once...
    lines = [_.strip().split() for _ in filename.open("r").readlines()]
    if not lines:
        msg = f"Empty delta file from nucmer, {filename}"
        raise ValueError(msg)
    for line in lines:
        if line[0] == "NUCMER":  # Skip headers
            continue
        # Lines starting with ">" indicate which sequences are aligned
        if line[0].startswith(">"):
            current_ref = line[0].strip(">")
            current_qry = line[1]
        # Lines with seven columns are alignment region headers:
        if len(line) == 7:  # noqa: PLR2004
            contains_alignment = True
            # Obtaining aligned regions needed to check for overlaps
            regions_ref[current_ref].append(
                tuple(sorted([int(line[0]), int(line[1])]))
            )  # aligned regions reference
            regions_qry[current_qry].append(
                tuple(sorted([int(line[2]), int(line[3])]))
            )  # aligned regions qry

            # Calculate aligned bases for each sequence
            ref_aln_lengths = abs(int(line[1]) - int(line[0])) + 1
            qry_aln_lengths = abs(int(line[3]) - int(line[2])) + 1
            aligned_bases += ref_aln_lengths + qry_aln_lengths

            # Calculate weighted identical bases
            sim_error += int(line[4])
            weighted_identical_bases += (ref_aln_lengths + qry_aln_lengths) - (
                2 * int(line[5])
            )

    # Calculate average %ID
    try:
        avrg_identity = float(weighted_identical_bases / aligned_bases)
    except ZeroDivisionError:
        avrg_identity = None  # Using Python's None to represent NULL

    # Return result based on whether alignments were found
    if not contains_alignment:
        return (None, None, None, None)
    return (
        get_aligned_bases_count(regions_qry),
        get_aligned_bases_count(regions_ref),
        avrg_identity,
        sim_error,
    )
