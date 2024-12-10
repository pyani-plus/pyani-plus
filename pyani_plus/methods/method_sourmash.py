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
"""Code to implement the sourmash Average Nucleotide Identity (ANI) method."""

# Set Up
import sys
from collections.abc import Iterator
from pathlib import Path

SCALED = 1000
KMER_SIZE = 31  # default


def parse_sourmash_compare_csv(
    compare_file: Path, filename_to_hash: dict[str, str]
) -> Iterator[tuple[str, str, float | None, float | None]]:
    """Parse sourmash compare all-vs-all CSV output.

    Assumes this is query-containment, and will infer the max-containment.

    Returns tuples of (query_hash, subject_hash, query-containment,
    max-containmenet) were the containment values of zero are mapped to
    None (to become null in the database).
    """
    query_containment: dict[tuple[str, str], float] = {}
    with compare_file.open() as handle:
        headers = handle.readline().rstrip("\n").split(",")
        if len(headers) != len(filename_to_hash):
            msg = (
                "Expected sourmash compare CSV to have"
                f" {len(filename_to_hash)} columns, not {len(headers)}"
            )
            sys.exit(msg)
        try:
            hashes = [filename_to_hash[Path(_).name] for _ in headers]
        except KeyError as err:
            msg = f"CSV file {compare_file} contained reference to {err!s} which is not in the run"
            sys.exit(msg)
        for row, query in enumerate(hashes):
            values = handle.readline().rstrip("\n").split(",")
            if values[row] != "1.0":
                # This could happen due to a bug in sourmash, or a glitch in the
                # parser if we're not looking at the matrix element we think we are?
                msg = (
                    f"Expected sourmash {query} vs self to be one, not {values[row]!r}"
                )
                raise ValueError(msg)
            for col, value in enumerate(values):
                query_containment[query, hashes[col]] = float(value)
    # Now that we have parsed the whole matrix,
    # can infer the max-containment
    for (query, subject), containment in query_containment.items():
        max_containment = max(containment, query_containment[subject, query])
        yield (
            query,
            subject,
            containment if containment else None,
            max_containment if max_containment else None,
        )


def parse_sourmash_manysearch_csv(
    manysearch_file: Path,
    filename_to_hash: dict[str, str],
    expected_pairs: set[tuple[str, str]],
) -> Iterator[tuple[str, str, float | None, float | None]]:
    """Parse sourmash-plugin-branchwater manysearch CSV output.

    Returns tuples of (query_hash, subject_hash, query-containment ANI
    estimate, max-containment ANI estimate).

    Any pairs not in the file are inferred to be failed alignments.
    """
    column_query = column_subject = column_query_cont = column_max_cont = 0
    with manysearch_file.open() as handle:
        line = handle.readline().rstrip("\n")
        headers = line.split(",")
        try:
            # In branchwater 0.9.11 column order varies between manysearch and pairwise
            column_query = headers.index("query_name")
            column_subject = headers.index("match_name")
            # Should be query_containment_ani, but that means transposing
            # the sourmash.csv compared to our identity matrix
            column_query_cont = headers.index("match_containment_ani")
            column_max_cont = headers.index("max_containment_ani")
        except ValueError:
            msg = f"ERROR - Missing expected fields in sourmash manysearch header: {line!r}"
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
                    f"Expected sourmash manysearch {filename_to_hash[values[column_query]]}"
                    f" vs self to be one, not {values[column_max_cont]!r}"
                )
                raise ValueError(msg)
            query_hash = filename_to_hash[values[column_query]]
            subject_hash = filename_to_hash[values[column_subject]]
            if (query_hash, subject_hash) in expected_pairs:
                expected_pairs.remove((query_hash, subject_hash))
            else:
                msg = f"Did not expect {query_hash} vs {subject_hash} in {manysearch_file.name}"
                raise ValueError(msg)
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
