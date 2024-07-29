import intervaltree
from collections import defaultdict
from pathlib import Path
from typing import Tuple


def parse_delta(filename, outputfile) -> Tuple[int, int, float, int]:
    """Return (reference alignment length, query alignment length, average identity, similarity erors)

    :param filename: Path to the input .delta file

    Calculates the aligned lengths for reference and query and average nucleotide
    identity, and returns the cumulative total for each as a tuple.

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

    current_ref, current_qry, raln_length, qaln_length, sim_error, avrg_ID = (
        None,
        None,
        0,
        0,
        0,
        0.0,
    )

    regions_ref = defaultdict(list)  # Hold a dictionary for query regions
    regions_qry = defaultdict(list)  # Hold a dictionary for query regions

    aligned_bases = []  # Hold a list for aligned bases for each sequence
    weighted_identical_bases = []  # Hold a list for weighted identical bases

    for line in [_.strip().split() for _ in Path(filename).open("r").readlines()]:
        if line[0] == "NUCMER":  # Skip headers
            continue
        # Lines starting with ">" indicate which sequences are aligned
        if line[0].startswith(">"):
            current_ref = line[0].strip(">")
            current_qry = line[1]
        # Lines with seven columns are alignment region headers:
        if len(line) == 7:
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
            aligned_bases.append(ref_aln_lengths)
            aligned_bases.append(qry_aln_lengths)

            # Calculate weighted identical bases
            sim_error += int(line[4])
            weighted_identical_bases.append(
                (ref_aln_lengths + qry_aln_lengths) - (2 * int(line[4]))
            )

    # Calculate average %ID
    avrg_ID = sum(weighted_identical_bases) / sum(aligned_bases)

    # Calculate total aligned bases (no overlaps)
    for seq_id in regions_qry:
        qry_tree = intervaltree.IntervalTree.from_tuples(regions_qry[seq_id])
        qry_tree.merge_overlaps(strict=False)
        for interval in qry_tree:
            qaln_length += interval.end - interval.begin + 1

    for seq_id in regions_ref:
        ref_tree = intervaltree.IntervalTree.from_tuples(regions_ref[seq_id])
        ref_tree.merge_overlaps(strict=False)
        for interval in ref_tree:
            raln_length += interval.end - interval.begin + 1
    values = (raln_length, qaln_length, avrg_ID, sim_error)

    with open(outputfile, "w") as file:
        for item in values:
            file.write(f"{item}\n")

    return (raln_length, qaln_length, avrg_ID, sim_error)


if __name__ == "__main__":
    import sys

    if sys.argv[1] == "parse_delta":
        input_file = sys.argv[2]
        output_file = sys.argv[3]
        parse_delta(input_file, output_file)
