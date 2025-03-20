# The MIT License
#
# Copyright (c) 2025 University of Strathclyde
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
"""Code to extract average nucleotide identity from an external alignment."""

import logging
from collections.abc import Callable, Iterator
from pathlib import Path

from pyani_plus import log_sys_exit

ASCII_GAP = ord("-")  # 45


def compute_external_alignment_column(  # noqa: PLR0913
    logger: logging.Logger,
    subject_hash: str,
    query_hashes: set[str],
    alignment: Path,
    mapping: Callable[[str], str],
    label: str,
) -> Iterator[tuple[str, str, float, int, int, float, float]]:
    """Parse a FASTA multiple alignment, and compute pairwise ANI values.

    Returns tuples of (query genome (hash), reference/subject genome (hash),
    ANI (float), alignment length (int), similarity errors (int), query
    coverage (float) and suject coverage (float).

    It will return two entries (a,b,...) and (b,a,...) for the off diagonals
    since this is a symmetric method.

    :param subject_hash: str, the genome MD5 hash to compute against
    :param query_hashes: Set[str], which query genome MD5 hashes to compare
    :param alignment: Path, path to the input FASTA alignment file
    :param mapping: Callable mapping alignment names to genome MD5 hashes
    :param label: str, description of the mapping for use in error messages
    """
    import numpy as np  # lazy import, although might be implicitly loaded already?

    from pyani_plus.utils import fasta_bytes_iterator

    # We could interpret the column number as the MSA ordering, but our internal
    # API by subject hash - and we can't assume the MSA is in any particular order.
    # Easiest way to solve this is two linear scans of the file, with a seek(0)
    with alignment.open("rb") as handle:
        subject_seq = subject_title = b""  # placeholder values
        s_non_gaps = subject_seq_gaps = -1  # placeholder values
        # Note loading the file in binary mode, will be working in bytes
        for query_title, query_seq in fasta_bytes_iterator(handle):
            query_hash = mapping(query_title.decode().split(None, 1)[0])
            if not query_hash:
                msg = (
                    f"Could not map {query_title.decode().split(None, 1)[0]} as {label}"
                )
                log_sys_exit(logger, msg)
            if query_hash == subject_hash:
                # for use in rest of the loop - an array of bytes!
                subject_seq = query_seq
                s_array = np.array(list(subject_seq), np.ubyte)
                s_non_gaps = s_array != ASCII_GAP
                subject_seq_gaps = subject_seq.count(b"-")
                subject_title = query_title  # for use in logging
                break
        else:
            msg = f"Did not find subject {subject_hash} in {alignment.name}"
            log_sys_exit(logger, msg)

        handle.seek(0)
        for query_title, query_seq in fasta_bytes_iterator(handle):
            query_hash = mapping(query_title.decode().split(None, 1)[0])
            if query_hash < subject_hash or query_hash not in query_hashes:
                # Exploiting symmetry to avoid double computation,
                # or not asked to compute this pairing (as already in the DB)
                continue
            if query_hash == subject_hash:
                # 100% identity and coverage, but need to calculate aln_length
                yield (
                    query_hash,
                    subject_hash,
                    1.0,
                    len(query_seq) - query_seq.count(b"-"),
                    0,
                    1.0,
                    1.0,
                )
            else:
                # Full calculation required
                if len(query_seq) != len(subject_seq):
                    msg = (
                        "Bad external-alignment, different lengths"
                        f" {len(query_seq)} and {len(subject_seq)}"
                        f" from {query_title.decode().split(None, 1)[0]}"
                        f" and {subject_title.decode().split(None, 1)[0]}``"
                    )
                    log_sys_exit(logger, msg)

                q_array = np.array(list(query_seq), np.ubyte)
                q_non_gaps = q_array != ASCII_GAP
                # & is AND
                # | is OR
                # ^ is XOR
                # ~ is NOT
                # e.g. ~(q_gaps | s_gaps) would be entries with no gaps
                naive_matches = q_array == s_array  # includes double gaps!
                matches = int((naive_matches & q_non_gaps).sum())
                one_gapped = q_non_gaps ^ s_non_gaps
                non_gap_mismatches = int((~naive_matches & ~one_gapped).sum())
                either_gapped = int(one_gapped.sum())
                del naive_matches, q_non_gaps, q_array

                # Now compute the alignment metrics from that
                query_cov = (matches + non_gap_mismatches) / (
                    len(query_seq) - query_seq.count(b"-")
                )
                subject_cov = (matches + non_gap_mismatches) / (
                    len(subject_seq) - subject_seq_gaps
                )
                aln_length = matches + non_gap_mismatches + either_gapped
                sim_errors = non_gap_mismatches + either_gapped

                yield (
                    query_hash,
                    subject_hash,
                    matches / aln_length,
                    aln_length,
                    sim_errors,
                    query_cov,
                    subject_cov,
                )
                # And the symmetric entry
                yield (
                    subject_hash,
                    query_hash,
                    matches / aln_length,
                    aln_length,
                    sim_errors,
                    subject_cov,
                    query_cov,
                )
