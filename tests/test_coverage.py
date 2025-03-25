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
"""Tests for coverage calculation.

These tests are intended to be run from the repository root using:

pytest -v
"""

import tempfile
from collections.abc import Callable
from pathlib import Path

from pyani_plus import db_orm, public_cli, setup_logger


def do_comparison(
    fasta_dir: Path,
    method: Callable,
    **kwargs: float | Path,
) -> db_orm.Run:
    """Execute an ANI method and return the run."""
    with tempfile.NamedTemporaryFile(suffix=".db") as tmp_db:
        method(
            database=tmp_db.name,
            fasta=fasta_dir,
            name="Artificial",
            create_db=True,
            **kwargs,
        )
        logger = setup_logger(None)
        session = db_orm.connect_to_db(logger, tmp_db.name)
        run = session.query(db_orm.Run).one()
        session.close()
        return run


def test_coverage(tmp_path: str) -> None:
    """Check comparison have expected identity and query-coverage."""
    tmp_dir = Path(tmp_path)
    seq_dir = tmp_dir / "fasta"
    seq_dir.mkdir()
    (seq_dir / "small.fasta").symlink_to(
        Path("tests/fixtures/MIBY01000005.fasta").resolve()
    )  # 154173fb8e7415ab45532a738572f957 - 7582 bp
    (seq_dir / "large.fasta").symlink_to(
        Path("tests/fixtures/MIBY01000011.fasta").resolve()
    )  # a0efc718e680e34d2f5c8f5d2286ca9c - 18001 bp
    with (seq_dir / "both.fasta").open("w") as handle:
        for file in (
            "tests/fixtures/MIBY01000005.fasta",
            "tests/fixtures/MIBY01000011.fasta",
        ):
            with Path(file).open() as in_handle:
                for line in in_handle:
                    handle.write(line)
        # Gives 7b6a6226ce00e52edca15565aa0d270d for both contigs
    assert len(list(seq_dir.glob("*.f*"))) == 3  # noqa: PLR2004

    checksums = (
        '["154173fb8e7415ab45532a738572f957",'  # small contig
        '"7b6a6226ce00e52edca15565aa0d270d",'  # both contigs
        '"a0efc718e680e34d2f5c8f5d2286ca9c"]'  # large contig
    )
    # Small vs large (and large vs small) should fail, giving null.
    # Self-vs-self ought to give 100% but we picked these two genomes
    # because they were odd for self-vs-self, so expect 99~100% identity:
    #
    # i.e. [[99~100%, 99~100%,    null],
    #       [99~100%, 99~100%, 99~100%],
    #       [   null, 99~100%, 99~100%]] for identity
    #
    # In terms of coverage, this should be 100% on the diagonal.
    # For the first row (small query) either 100% (vs self of the both file)
    # or null (vs the larger contig). For the middle row (both contigs),
    # coverage will be ~30% (vs small), 100% (vs self), and ~70% (vs large).
    # Finally, for the last row (large query), coverage will be null (vs
    # small), 100% (vs both), and 100% (vs self).
    #
    # i.e. [[100%, 100%, null],
    #       [~30%, 100%, ~70%],
    #       [null, 100%, 100%]] for query_cov
    #
    # Using the contig lengths, small is 7582 bp, large is 18001 bp, gives
    # the expected coverages values vs the combined file of 29.6% and 70.4%.
    #
    # Comparing the JSON serialised dataframe directly from the DB

    run = do_comparison(seq_dir, public_cli.cli_anim)
    assert run.df_identity == (
        "{"
        f'"columns":{checksums},"index":{checksums},"data":'
        "[[0.996307043,0.996307043,null],[0.996307043,0.9989055232,1.0],[null,1.0,1.0]]}"
    )  # all 99~100% except two null values
    assert run.df_cov_query == (
        "{"
        f'"columns":{checksums},"index":{checksums},"data":'
        "[[1.0,1.0,null],[0.2963686823,1.0,0.7036313177],[null,1.0,1.0]]}"
    )  # expected pattern of 100%, 30%, 70% or null.

    run = do_comparison(seq_dir, public_cli.cli_dnadiff)
    assert run.df_identity == (
        "{"
        f'"columns":{checksums},"index":{checksums},"data":'
        "[[1.0,1.0,null],[1.0,1.0,1.0],[null,1.0,1.0]]}"
    )  # all 100% except two null values
    assert run.df_cov_query == (
        "{"
        f'"columns":{checksums},"index":{checksums},"data":'
        "[[1.0,1.0,null],[0.2963686823,1.0,0.7036313177],[null,1.0,1.0]]}"
    )

    run = do_comparison(seq_dir, public_cli.cli_anib)
    assert run.df_identity == (
        "{"
        f'"columns":{checksums},"index":{checksums},"data":'
        "[[1.0,1.0,null],[1.0,1.0,1.0],[null,1.0,1.0]]}"
    )  # all 100% except two null values
    assert run.df_cov_query == (
        "{"
        f'"columns":{checksums},"index":{checksums},"data":'
        "[[1.0,1.0,null],[0.2963686823,1.0,0.7036313177],[null,1.0,1.0]]}"
    )

    # Deliberately trying some non-default settings with fastANI
    run = do_comparison(
        seq_dir,
        public_cli.cli_fastani,
        kmersize=15,
        fragsize=2000,
        minmatch=0.15,
    )
    # The large contig is a non-100% self test case for fastANI
    assert run.df_identity == (
        "{"
        f'"columns":{checksums},"index":{checksums},"data":'
        "[[1.0,1.0,null],[1.0,0.99997,0.999959],[null,0.999959,0.999959]]}"
    )  # all 99~100% except two null values
    assert run.df_cov_query == (
        "{"
        f'"columns":{checksums},"index":{checksums},"data":'
        "[[1.0,1.0,null],[0.25,1.0,0.75],[null,1.0,1.0]]}"
    )  # 25% and 75% rather than 30% and 70% expected via bp

    # Doesn't "work" with default scaling - nulls except for diagonal,
    run = do_comparison(seq_dir, public_cli.cli_sourmash, scaled=50, cache=tmp_dir)
    assert run.df_identity == (
        "{"
        f'"columns":{checksums},"index":{checksums},"data":'
        "[[1.0,1.0,null],[1.0,1.0,1.0],[null,1.0,1.0]]}"
    )  # all 100% except two null values
    # get k-mer query-containment 96.2% and 98.8% where query-coverage was 30% and 70%
    assert run.df_cov_query == (
        "{"
        f'"columns":{checksums},"index":{checksums},"data":'
        "[[1.0,1.0,null],[0.9622440235,1.0,0.9884105907],[null,1.0,1.0]]}"
    )
