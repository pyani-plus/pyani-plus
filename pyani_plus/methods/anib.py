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

import gzip
from pathlib import Path

from pyani_plus.utils import fasta_bytes_iterator

FRAGSIZE = 1020  # Default ANIb fragment size
MIN_COVERAGE = 0.7
MIN_IDENTITY = 0.3

# std = qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
# We do NOT use the standard 12 columns, nor a custom 15 cols as per old pyANI,
# but a minimal set of only 7 (the 6 used fields plus subject id for debugging).
BLAST_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "qlen", "gaps"]

# Precompute the column indexes once
BLAST_COL_QUERY = BLAST_COLUMNS.index("qseqid")  # alternative to qaccver in std
BLAST_COL_PIDENT = BLAST_COLUMNS.index("pident")  # in std 12 columns
BLAST_COL_LENGTH = BLAST_COLUMNS.index("length")  # in std 12 columns
BLAST_COL_MISMATCH = BLAST_COLUMNS.index("mismatch")  # in std 12 columns
BLAST_COL_QLEN = BLAST_COLUMNS.index("qlen")
BLAST_COL_GAPS = BLAST_COLUMNS.index("gaps")


def fragment_fasta_file(
    filename: Path, fragmented_fasta: Path, fragsize: int = FRAGSIZE
) -> None:
    """Fragment FASTA file into subsequences of up to the given size.

    Any remainder is taken as is (1 <= length < fragsize).

    Accepts gzipped files as input.
    """
    with (
        (
            gzip.open(filename, "rb")
            if filename.suffix == ".gz"
            else filename.open("rb")
        ) as in_handle,
        fragmented_fasta.open("wb") as out_handle,
    ):
        count = 0
        for title, seq in fasta_bytes_iterator(in_handle):
            index = 0
            while index < len(seq):
                count += 1
                fragment = seq[index : index + fragsize]
                out_handle.write(b">frag%05i %s\n" % (count, title))
                # Now line wrap at 60 chars
                for i in range(0, len(fragment), 60):
                    out_handle.write(fragment[i : i + 60] + b"\n")
                index += fragsize
    if not count:
        msg = f"No sequences found in {filename}"
        raise ValueError(msg)


def parse_blastn_file(blastn: Path) -> tuple[float | None, int | None, int | None]:
    """Extract the ANI etc from a blastn output file using the ANIb method.

    Parses the BLAST tabular output file, taking only rows with a query coverage
    over 70% and percentage identity over 30%, and deduplicating.

    Returns mean percentage identity of all BLAST alignments passing thresholds
    (treated as zero rather than NaN if there are no accepted alignments), the
    total alignment length, and total similarity errors (mismatches and gaps).

    >>> fname = "tests/fixtures/viral_example/intermediates/ANIb/MGV-GENOME-0264574_vs_MGV-GENOME-0266457.tsv"
    >>> identity, length, sim_errors = parse_blastn_file(Path(fname))
    >>> print(
    ...     f"Identity {100 * identity:0.1f}% over length {length} with {sim_errors} errors"
    ... )
    Identity 99.5% over length 39169 with 215 errors

    We expect 100% identity for a self comparison (but this is not always true):

    >>> fname = "tests/fixtures/viral_example/intermediates/ANIb/MGV-GENOME-0264574_vs_MGV-GENOME-0264574.tsv"
    >>> identity, length, sim_errors = parse_blastn_file(Path(fname))
    >>> print(
    ...     f"Identity {100 * identity:0.1f}% over length {length} with {sim_errors} errors"
    ... )
    Identity 100.0% over length 39253 with 0 errors

    Returns (None, None, None) for no acceptable alignment.
    """
    total_pid_100 = 0.0
    total_count = 0
    total_aln_length = 0
    total_sim_errors = 0

    prev_query = ""
    with blastn.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(BLAST_COLUMNS):
                msg = (
                    f"Found {len(fields)} columns in {blastn}, not {len(BLAST_COLUMNS)}"
                )
                raise ValueError(msg)
            query = fields[BLAST_COL_QUERY]
            if not query.startswith("frag"):
                msg = f"BLAST output should be using fragmented queries, not {query}"
                raise ValueError(msg)
            blast_gaps = int(fields[BLAST_COL_GAPS])
            ani_alnlen = int(fields[BLAST_COL_LENGTH]) - blast_gaps
            blast_mismatch = int(fields[BLAST_COL_MISMATCH])
            ani_query_coverage = ani_alnlen / int(fields[BLAST_COL_QLEN])
            # Can't use float(values["pident"])/100, this is relative to alignment length
            ani_pid = (ani_alnlen - blast_mismatch) / int(fields[BLAST_COL_QLEN])

            # Now apply filters - should these be parameters?
            # And if there are multiple hits for this query, take first (best) one
            if (
                ani_query_coverage > MIN_COVERAGE
                and ani_pid > MIN_IDENTITY
                and prev_query != query
            ):
                total_aln_length += ani_alnlen
                total_sim_errors += blast_mismatch + blast_gaps
                # Not using ani_pid but BLAST's pident - see note below:
                total_pid_100 += float(fields[BLAST_COL_PIDENT])
                total_count += 1
                prev_query = query  # to detect multiple hits for a query
    # NOTE: Could warn about empty BLAST file using if prev_query is None:

    # NOTE: We report the mean of blastn's pident for concordance with JSpecies
    # Despite this, the concordance is not exact. Manual inspection during
    # the original pyANI development indicated that a handful of fragments
    # are differentially filtered out in JSpecies and here. This is often
    # on the basis of rounding differences (e.g. coverage being close to 70%).
    return (
        total_pid_100 / (total_count * 100) if total_count else None,
        total_aln_length if total_count else None,
        total_sim_errors if total_count else None,
    )
