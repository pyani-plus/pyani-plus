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
MIN_COVERAGE = 0.7
MIN_IDENTITY = 0.3

# We do NOT use the standard 12 columns, but a custom 15 cols as per old pyANI
# std = qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
BLAST_COLUMNS = (
    "qseqid sseqid length mismatch pident nident"
    " qlen slen qstart qend sstart send positive ppos gaps"
).split()


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
        if not count:
            msg = f"No sequences found in {filename}"
            raise ValueError(msg)
        fragmented_files.append(frag_filename)
    return fragmented_files


def parse_blastn_file(blastn: Path) -> tuple[float, int, int]:
    """Extract the ANI etc from a blastn output file using the ANIb method.

    Parses the BLAST tabular output file, taking only rows with a query coverage
    over 70% and percentage identity over 30%, and deduplicating.

    Returns mean percentage identity of all BLAST alignments passing thresholds
    (treated as zero rather than NaN if there are no accepted alignments), the
    total alignment length, and total similarity errors (mismatches and gaps).

    >>> fname = (
    ...     "tests/fixtures/anib/blastn/MGV-GENOME-0264574_vs_MGV-GENOME-0266457.tsv"
    ... )
    >>> identity, length, sim_errors = parse_blastn_file(Path(fname))
    >>> print(
    ...     f"Identity {100*identity:0.1f}% over length {length} with {sim_errors} errors"
    ... )
    Identity 99.5% over length 39169 with 215 errors

    We expect 100% identity for a self comparison (but this is not always true):

    >>> fname = (
    ...     "tests/fixtures/anib/blastn/MGV-GENOME-0264574_vs_MGV-GENOME-0264574.tsv"
    ... )
    >>> identity, length, sim_errors = parse_blastn_file(Path(fname))
    >>> print(
    ...     f"Identity {100*identity:0.1f}% over length {length} with {sim_errors} errors"
    ... )
    Identity 100.0% over length 39253 with 0 errors
    """
    total_aln_length = 0
    total_sim_errors = 0
    all_pid: list[float] = []

    prev_query = ""
    with blastn.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(BLAST_COLUMNS):
                msg = (
                    f"Found {len(fields)} columns in {blastn}, not {len(BLAST_COLUMNS)}"
                )
                raise ValueError(msg)
            if not fields[0].startswith("frag"):
                msg = (
                    f"BLAST output should be using fragmented queries, not {fields[0]}"
                )
                raise ValueError(msg)
            values = dict(zip(BLAST_COLUMNS, fields, strict=False))
            blast_alnlen = int(values["length"])
            blast_gaps = int(values["gaps"])
            ani_alnlen = blast_alnlen - blast_gaps
            blast_mismatch = int(values["mismatch"])
            ani_alnids = ani_alnlen - blast_mismatch
            ani_query_coverage = ani_alnlen / int(values["qlen"])
            # Can't use float(values["pident"])/100, this is relative to alignment length
            ani_pid = ani_alnids / int(values["qlen"])

            # Now apply filters - should these be parameters?
            # And if there are multiple hits for this query, take first (best) one
            if (
                ani_query_coverage > MIN_COVERAGE
                and ani_pid > MIN_IDENTITY
                and prev_query != fields[0]
            ):
                total_aln_length += ani_alnlen
                total_sim_errors += blast_mismatch + blast_gaps
                # Not using ani_pid but BLAST's pident - see note below:
                all_pid.append(float(values["pident"]) / 100)
                prev_query = fields[0]  # to detect multiple hits for a query
    # NOTE: Could warn about empty BLAST file using if prev_query is None:

    # NOTE: We report the mean of blastn's pident for concordance with JSpecies
    # Despite this, the concordance is not exact. Manual inspection during
    # the original pyANI development indicated that a handful of fragments
    # are differentially filtered out in JSpecies and here. This is often
    # on the basis of rounding differences (e.g. coverage being close to 70%).
    return (
        sum(all_pid) / len(all_pid) if all_pid else 0,
        total_aln_length,
        total_sim_errors,
    )
