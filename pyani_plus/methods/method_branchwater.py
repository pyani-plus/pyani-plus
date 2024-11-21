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
"""Code to implement the sourmash-plugin-branchwater ANI method."""

# Set Up
import sys
from collections.abc import Iterator
from pathlib import Path


def parse_sourmash_pairwise_csv(
    pairwise_file: Path, filename_to_hash: dict[str, str]
) -> Iterator[tuple[str, str, float]]:
    """Parse sourmash-plugin-branchwater pairwise CSV output.

    Returns tuples of (query_hash, subject_hash, estimated ANI).
    """
    column_query = 0
    column_subject = 2
    column_ani = 14  # use mode, here assuming max-containment!
    with pairwise_file.open() as handle:
        headers = handle.readline().rstrip("\n").split(",")
        if headers != [
            "query_name",  # column 0
            "query_md5",  # this is the sig file's MD5
            "match_name",  # column 2
            "match_md5",  # this is the sig file's MD5
            "containment",
            "max_containment",
            "jaccard",
            "intersect_hashes",
            "ksize",
            "scaled",
            "moltype",
            "query_containment_ani",
            "match_containment_ani",
            "average_containment_ani",
            "max_containment_ani",  # column 14
        ]:
            msg = (
                "ERROR: Wrong header for sourmash-plugin-breakwater pairwise CSV output"
            )
            sys.exit(msg)
        for line in handle:
            values = line.rstrip("\n").split(",")
            if len(values) != 15:  # noqa: PLR2004
                msg = f"sourmash manysearch line did not have 15 fields: {line!r}"
                raise ValueError(msg)
            yield (
                filename_to_hash[values[column_query]],
                filename_to_hash[values[column_subject]],
                float(values[column_ani]),
            )


def parse_sourmash_manysearch_csv(
    manysearch_file: Path, filename_to_hash: dict[str, str]
) -> Iterator[tuple[str, str, float]]:
    """Parse sourmash-plugin-branchwater manysearch CSV output.

    Returns tuples of (query_hash, subject_hash, estimated ANI).
    """
    column_query = 0
    column_subject = 2
    column_ani = 14  # use mode, here assuming max-containment!
    with manysearch_file.open() as handle:
        line = handle.readline().rstrip("\n")
        headers = line.split(",")
        if headers != [
            "query_name",  # column 0
            "query_md5",  # this is the hash of the sig
            "match_name",  # column 2
            "containment",
            "intersect_hashes",
            "ksize",
            "scaled",
            "moltype",
            "match_md5",
            "jaccard",
            "max_containment",
            "query_containment_ani",
            "match_containment_ani",
            "average_containment_ani",
            "max_containment_ani",  # column 14
        ]:
            msg = f"ERROR: Wrong header for sourmash-plugin-breakwater pairwise CSV output:\n{line!r}"
            sys.exit(msg)
        for line in handle:
            values = line.rstrip("\n").split(",")
            if not values:
                continue
            if len(values) != 15:  # noqa: PLR2004
                msg = f"sourmash manysearch line did not have 15 fields: {line!r}"
                raise ValueError(msg)
            if (
                values[column_query] == values[column_subject]
                and values[column_ani] != "1.0"
            ):
                msg = f"sourmash self-vs-self was {values[column_ani]!r} not one"
                raise ValueError(msg)
            yield (
                filename_to_hash[values[column_query]],
                filename_to_hash[values[column_subject]],
                float(values[column_ani]),
            )
