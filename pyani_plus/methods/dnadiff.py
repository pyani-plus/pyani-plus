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
"""Code to implement the dnadiff Average Nucleotide Identity (ANI) method.

This script calculates ANI using the method implemented by the dnadiff
software, available at https://github.com/garviz/MUMmer/blob/master/dnadiff.

All input FASTA format files are compared against each other in both
directions using NUCmer with the --maxmatch parameter. To replicate the
AlignedBases value reported by dnadiff software, delta-filter (-m parameters),
show-diff (-qH parameters), and show-coords (-rclTH parameters) are also run.
All outputs will be stored in a specified output directory.

Outputs generated by show-coords (.mcoords) and show-diff (.qdiff) files are
parsed to obtain values for the total number of aligned bases (AlignedBases),
average nucletide identity (ANI), and genome coverage for each pairwise
comparison.

To report the total number of aligned bases, we start by parsing .mcoords file,
which has 13 columns, with each line representing a single alignment. See
https://mummer.sourceforge.net/manual/#coords for specification:

```
85	37713	1	37636	37629	37636	99.43	39253	39594	95.86	95.05	MGV_MGV-GENOME-0264574	MGV_MGV-GENOME-0266457
17709	39253	17626	39176	21545	21551	99.97	39253	39594	54.89	54.43	MGV_MGV-GENOME-0264574	MGV_MGV-GENOME-0266457
```

- col_0:  start of the alignment region in the reference sequence
- col_1:  end of the alignment region in the reference sequence
- col_2: start of the alignment region in the query sequence
- col_3: end of the alignment region in the query sequence
- col_4: length of the alignment region in the reference sequence
- col_5: length of the alignment region in the query sequence
- col_6: percent identity of the alignment
- col_7: length of the reference sequence
- col_8: length of the query sequence
- col_9: percent alignment coverage in the reference sequence
- col_10:  percent alignment coverage in the query sequence
- col_11: reference sequence ID
- col_12: query sequence ID

We check if the sequence from the query genome (col_12) has been seen before.
If it has not, we update the total number of aligned bases by adding the total
length of that sequence (col_8). At this point, the final value will be the
total length of the alignment, including gaps (eg. 39594 for the query sequence).

In the next step, we look at the .qdiff files generated by show-coords, which
output a list of structural differences for each sequence in the reference
and query, sorted by position.

```
MGV_MGV-GENOME-0266457	GAP	37637	17625	-20011	-20005	-6
MGV_MGV-GENOME-0266457	BRK	39177	39594	418
```
IDQ GAP gap-start gap-end gap-length-Q gap-length-R gap-diff
IDQ DUP dup-start dup-end dup-length
IDQ BRK gap-start gap-end gap-length
IDQ JMP gap-start gap-end gap-length
IDQ INV gap-start gap-end gap-length
IDQ SEQ gap-start gap-end gap-length prev-sequence next-sequence

The total number of gaps in the alignments is calculated by considering
only structural differences other than duplications (DUP) and the positive gap
values (gap > 0) in the col_4 (NOTE: column numbers are from 0 to 7). The final
number of non-redundant bases in the query genome aligned to the reference genome
is calculated by subtracting the lengths of the gaps from the previously
calculated total length of the alignment, including gaps (eg. 39594 - 418 = 39176).

We report genome coverage as the percentage of the query genome that is aligned
to the reference genome. This is calculated as the number of query aligned
bases without gaps, divided by the total length of the query genome
(eg. (39594 / 39176) = 98.94%).

The average nucleotide identity is reported from the .mcoords file and is
calculated as the total sum of identities divided by the total sum ofalignment
lengths, multiplied by 100. The sum of identity for a single alignment is
calculated as follows:
(percent identity of the alignment [col_6] / 100) * (QRY + REF alignment length (col_4 + col_5)),

where the alignment sum for a single alignments is calculated as:

QRY + REF alignment length (col_4 + col_5).
"""  # noqa: E501

# Set Up
from pathlib import Path

import pandas as pd


def parse_mcoords(mcoords_file: Path) -> tuple[float | None, int | None]:
    """Parse mcoords file and return avg ID% and number of aligned bases with gaps (QRY).

    :parama mcoords_file: Path to the mcoords_file
    """
    mcoords = pd.read_csv(
        Path(mcoords_file),
        sep="\t",
        names=[f"col_{_}" for _ in range(13)],
    )
    if mcoords.empty:
        return (None, None)
    sum_identity = 0.0  # should be an int, but we can only estimate it
    sum_alignment_lengths = 0
    seen_ref_seq = set()
    aligned_bases_with_gaps = 0

    # Calculate average nucleotide identity
    for _index, row in mcoords.iterrows():
        row_length = int(row["col_4"]) + int(row["col_5"])
        sum_identity += float(row["col_6"]) * row_length / 100
        sum_alignment_lengths += row_length

        # Get number of aligned bases with gaps
        if row["col_12"] not in seen_ref_seq:
            aligned_bases_with_gaps += int(row["col_8"])
            seen_ref_seq.add(row["col_12"])
    return (sum_identity / sum_alignment_lengths, aligned_bases_with_gaps)


def parse_qdiff(qdiff_file: Path) -> int | None:
    """Parse qdiff and return number of gaps in the QRY alignment.

    :param qdiff_file: Path to the qdiff file
    """
    qdiff = pd.read_csv(
        Path(qdiff_file),
        sep="\t",
        names=[f"col_{_}" for _ in range(7)],
    )
    if qdiff.empty:
        return None
    gap_lengths = 0
    for _index, row in qdiff.iterrows():
        gap = row["col_4"]
        if row["col_1"] != "DUP" and gap > 0:
            gap_lengths += gap

    return gap_lengths
