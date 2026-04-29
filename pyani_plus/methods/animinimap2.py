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
"""Code to implement ANI but with minimap2 rather than MUMmer.

The method is inspired by the ANIm method described in Richter et al (2009),
Proc Natl Acad Sci USA 106: 19126-19131 doi:10.1073/pnas.0906412106.

All input FASTA format files are compared against each other, pairwise, using
Heng Li's minimap2:

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
Bioinformatics, 34:3094-3100. https://doi.org/10.1093/bioinformatics/bty191

Li, H. (2021). New strategies to improve minimap2 alignment accuracy.
Bioinformatics, 37:4572-4574. https://doi.org/10.1093/bioinformatics/btab705

For efficiency, each reference genome is first indexed once, and the index reused.
The PAF output is parsed to obtain an alignment length and exact match count.
These are processed to give matrices of aligned sequence lengths, average
nucleotide identity (ANI) percentages, and minimum aligned percentage (of whole
genome) for each pairwise comparison.
"""

from collections import defaultdict
from pathlib import Path

from pyani_plus.methods.anim import get_aligned_bases_count
from pyani_plus.public_cli_args import EnumPresetMinimap2

DEFAULT_PRESET = EnumPresetMinimap2.asm10

PAF_COL_QUERY_NAME = 0
PAF_COL_QUERY_LENGTH = 1
PAF_COL_QUERY_START = 2
PAF_COL_QUERY_END = 3
PAF_COL_RELATIVE_STRAND = 4
PAF_COL_TARGET_NAME = 5
PAF_COL_TARGET_LENGTH = 6
PAF_COL_TARGET_START = 7
PAF_COL_TARGET_END = 8
PAF_COL_MATCHES = 9
PAF_COL_ALN_BLOCK_LEN = 10
PAF_COL_MAPPING_QUALITY = 11


def parse_minimap2_paf_file(filename: Path) -> tuple[int, int, float | None]:
    """Return (reference alignment length, query alignment length, average identity).

    :param filename: Path to a minimap2 PAF output file.

    Calculates similarity errors and the aligned lengths for reference
    and query and average nucleotide identity, and returns the cumulative
    total for each as a tuple.

    The minimap PAF file format contains 12 numbers in each line, some of which
    are of interest: https://github.com/lh3/miniasm/blob/master/PAF.md

    We report ANI identity by finding an average across all alignments using
    the following formula:

    sum of identical bases * 2 / sum of aligned bases from each fragment

    To calculate alignment lengths, we extract the regions of each alignment
    (separately for query or reference) provided in the PAF file and merge
    the overlapping regions with IntervalTree. Then, we calculate the total
    sum of all aligned regions.
    """
    regions_ref = defaultdict(list)  # Hold a dictionary for query regions
    regions_qry = defaultdict(list)  # Hold a dictionary for query regions

    aligned_bases = 0  # Hold a count of aligned bases for each sequence
    identical_bases = 0  # Hold a count of identical bases

    # Ideally we wouldn't read the whole file into memory at once...
    lines = [_.strip().split() for _ in filename.open("r").readlines()]
    if not lines:
        msg = f"Empty PAF file from minimap2, {filename}"
        raise ValueError(msg)
    for line in lines:
        if line[PAF_COL_MAPPING_QUALITY] == "0":
            # Skip failed alignments - see also the SAM output
            continue
        current_ref = line[PAF_COL_TARGET_NAME]
        current_qry = line[PAF_COL_QUERY_NAME]

        # Obtaining aligned regions needed to check for overlaps
        regions_ref[current_ref].append(
            tuple(
                sorted([int(line[PAF_COL_TARGET_START]), int(line[PAF_COL_TARGET_END])])
            )
        )  # aligned regions reference
        regions_qry[current_qry].append(
            tuple(
                sorted([int(line[PAF_COL_QUERY_START]), int(line[PAF_COL_QUERY_END])])
            )
        )  # aligned regions query

        # Calculate aligned bases for each sequence
        ref_aln_lengths = (
            abs(int(line[PAF_COL_TARGET_END]) - int(line[PAF_COL_TARGET_START])) + 1
        )
        qry_aln_lengths = (
            abs(int(line[PAF_COL_QUERY_END]) - int(line[PAF_COL_QUERY_START])) + 1
        )
        aligned_bases += ref_aln_lengths + qry_aln_lengths

        # Calculate weighted identical bases
        identical_bases += int(line[PAF_COL_MATCHES])

    # Calculate average %ID
    try:
        avrg_identity = float(identical_bases * 2 / aligned_bases)
    except ZeroDivisionError:
        avrg_identity = None  # Using Python's None to represent NULL

    return (
        get_aligned_bases_count(regions_qry),
        get_aligned_bases_count(regions_ref),
        avrg_identity,
    )
