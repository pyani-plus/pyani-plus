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
"""Code to wrap the fastANI average nucleotide identity method."""

from collections.abc import Iterator
from pathlib import Path

KMER_SIZE = 16  # Default fastANI k-mer size (max 16)
FRAG_LEN = 3000  # Default fastANI fragment length, mapped to our --fragsize
MIN_FRACTION = 0.2  # Default fastANI min fraction
MAX_RATIO_DIFF = 10.0  # Default fastANI maximum ratio difference


def parse_fastani_file(
    filename: Path,
    filename_to_hash: dict[str, str],
    expected_pairs: set[tuple[str, str]],
) -> Iterator[tuple[str, str, float | None, int | None, int | None]]:
    """Parse a multi-line fastANI output file extracting key fields as a tuple.

    Returns tuples of (query genome (hash), reference/subject genome (hash),
    ANI estimate (float or None), orthologous matches (int or None), and
    sequence fragments (float or None)).

    :param filename: Path, path to the input file
    :param filename_to_hash: Dict, mapping filename stems to MD5 hashes
    :param expected_pairs: Set, expected (query hash, subject hash) pairs

    Extracts the ANI estimate (which we return in the range 0 to 1), the
    number of orthologous matches (int), and the number of sequence
    fragments considered from the fastANI output file (int).

    Failed alignments are inferred by not appearing in the file), and get
    None values.

    Example fastANI comparison, three queries vs one reference:

    >>> mapping = {
    ...     "MGV-GENOME-0264574.fas": "A",
    ...     "MGV-GENOME-0266457.fna": "B",
    ...     "OP073605.fasta": "C",
    ... }
    >>> fname = Path(
    ...     "tests/fixtures/viral_example/intermediates/fastANI/all_vs_OP073605.fastani"
    ... )
    >>> expected = {(_, "C") for _ in ("A", "B", "C")}
    >>> for query, subject, ani, matches, frags in parse_fastani_file(
    ...     fname, mapping, expected
    ... ):
    ...     print(f"{query} vs {subject} gave {100 * ani:0.1f}%")
    A vs C gave 99.9%
    B vs C gave 99.5%
    C vs C gave 100.0%

    All the query-subject pairs in the file must be in expected_pairs:

    >>> for query, subject, ani, matches, frags in parse_fastani_file(
    ...     fname, mapping, {("A", "B"), ("B", "A")}
    ... ):
    ...     print(
    ...         f"{query} vs {subject} gave {100 * ani:0.1f}%"
    ...         if ani
    ...         else f"{query} vs {subject} gave {ani}"
    ...     )
    Traceback (most recent call last):
      ...
    ValueError: Did not expect A vs C in all_vs_OP073605.fastani

    Any pairs not in the file are inferred to be failed alignments, even
    for empty files:

    >>> fname = Path("/dev/null")
    >>> for query, subject, ani, matches, frags in parse_fastani_file(
    ...     fname, mapping, {("A", "B")}
    ... ):
    ...     print(query, subject, ani, matches, frags)
    A B None None None
    """
    with filename.open() as handle:
        for line in handle:
            parts = line.strip().split("\t")
            query_hash = filename_to_hash[Path(parts[0]).name]
            subject_hash = filename_to_hash[Path(parts[1]).name]
            if (query_hash, subject_hash) in expected_pairs:
                expected_pairs.remove((query_hash, subject_hash))
            else:
                msg = (
                    f"Did not expect {query_hash} vs {subject_hash} in {filename.name}"
                )
                raise ValueError(msg)
            yield (
                query_hash,
                subject_hash,
                0.01 * float(parts[2]),
                int(parts[3]),
                int(parts[4]),
            )
    # Even if the file was empty, we infer any remaining pairs
    # are failed alignments:
    for query_hash, subject_hash in expected_pairs:
        yield query_hash, subject_hash, None, None, None
