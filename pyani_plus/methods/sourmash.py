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
"""Code to implement the sourmash Average Nucleotide Identity (ANI) method."""

# Set Up
import sys
from collections.abc import Iterator
from pathlib import Path

SCALED = 1000
KMER_SIZE = 31  # default


def parse_sourmash_manysearch_csv(
    manysearch_file: Path,
    expected_pairs: set[tuple[str, str]],
) -> Iterator[tuple[str, str, float | None, float | None]]:
    """Parse sourmash-plugin-branchwater manysearch CSV output.

    Returns tuples of (query_hash, subject_hash, query-containment ANI
    estimate, max-containment ANI estimate).

    Any pairs not in the file are inferred to be failed alignments.  As we are
    using reporting threshold zero, this only applies to corner cases where
    there are no common k-mers.
    """
    column_query = column_subject = column_query_cont = column_max_cont = 0
    with manysearch_file.open() as handle:
        line = handle.readline().rstrip("\n")
        headers = line.split(",")
        try:
            # In branchwater 0.9.11 column order varies between manysearch and pairwise
            column_query = headers.index("query_name")
            column_subject = headers.index("match_name")
            column_query_cont = headers.index("query_containment_ani")
            column_max_cont = headers.index("max_containment_ani")
        except ValueError:
            msg = f"ERROR - Missing expected fields in sourmash manysearch header, found: {line!r}"
            sys.exit(msg)
        for line in handle:
            line = line.rstrip("\n")  # noqa: PLW2901
            if not line:
                continue
            values = line.split(",")
            if (
                values[column_query] == values[column_subject]
                and values[column_max_cont] != "1.0"
            ):
                msg = (
                    f"Expected sourmash manysearch {values[column_query]}"
                    f" vs self to be one, not {values[column_max_cont]!r}"
                )
                raise ValueError(msg)
            query_hash = values[column_query]
            subject_hash = values[column_subject]
            if (query_hash, subject_hash) in expected_pairs:
                expected_pairs.remove((query_hash, subject_hash))
            else:
                msg = f"Did not expect {query_hash} vs {subject_hash} in {manysearch_file.name}"
                raise ValueError(msg) from None
            yield (
                query_hash,
                subject_hash,
                float(values[column_query_cont]),
                float(values[column_max_cont]),
            )
    # Even if the file was empty (bar the header),
    # we infer any remaining pairs are failed alignments:
    for query_hash, subject_hash in expected_pairs:
        yield query_hash, subject_hash, None, None
