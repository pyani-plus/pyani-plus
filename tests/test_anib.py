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
"""Tests for the ANIb implementation.

These tests are intended to be run from the repository root using:

pytest -v
"""

import filecmp
from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, setup_logger, tools
from pyani_plus.methods import anib

from . import get_matrix_entry


def test_bad_path(tmp_path: str) -> None:
    """Confirm giving an empty path etc fails."""
    with pytest.raises(
        FileNotFoundError, match="No such file or directory: '/does/not/exist'"
    ):
        anib.parse_blastn_file(Path("/does/not/exist"))

    with pytest.raises(
        FileNotFoundError, match="No such file or directory: '/does/not/exist'"
    ):
        anib.fragment_fasta_file(
            Path("/does/not/exist"), Path(tmp_path) / "frags.fna", 1020
        )


def test_empty_path(tmp_path: str) -> None:
    """Confirm fragmenting an empty path fails."""
    with pytest.raises(ValueError, match="No sequences found in /dev/null"):
        anib.fragment_fasta_file(Path("/dev/null"), Path(tmp_path) / "frags.fna", 1020)


def test_parse_blastn_empty() -> None:
    """Check parsing of empty BLASTN tabular file."""
    assert anib.parse_blastn_file(Path("/dev/null")) == (None, None, None)


def test_parse_blastn_bad(input_genomes_tiny: Path) -> None:
    """Check parsing something which isn't a BLAST TSV files fails."""
    with pytest.raises(
        ValueError, match="Found 1 columns in .*/MGV-GENOME-0264574.fas, not 7"
    ):
        anib.parse_blastn_file(input_genomes_tiny / "MGV-GENOME-0264574.fas")


def test_parse_blastn_bad_query(tmp_path: str) -> None:
    """Check parsing TSV not using frag#### fails."""
    fake = Path(tmp_path) / "fake.tsv"
    with fake.open("w") as handle:
        handle.write("\t".join(["bad-query", "subject"] + ["0"] * 5))
    with pytest.raises(
        ValueError,
        match="BLAST output should be using fragmented queries, not bad-query",
    ):
        anib.parse_blastn_file(fake)


def test_parse_blastn(input_genomes_tiny: Path) -> None:
    """Check parsing of BLASTN tabular file."""
    # Function returns tuple of mean percentage identity, total alignment length, and
    # total mismatches/gaps:
    assert anib.parse_blastn_file(
        input_genomes_tiny
        / "intermediates/ANIb/MGV-GENOME-0264574_vs_MGV-GENOME-0266457.tsv"
    ) == (0.9945938461538462, 39169, 215)
    # $ md5sum tests/fixtures/viral_example/*.f*
    # 689d3fd6881db36b5e08329cf23cecdd tests/fixtures/viral_example/MGV-GENOME-0264574.fas
    # 78975d5144a1cd12e98898d573cf6536 tests/fixtures/viral_example/MGV-GENOME-0266457.fna
    # 5584c7029328dc48d33f95f0a78f7e57 tests/fixtures/viral_example/OP073605.fasta
    #
    # Example is query 689d3fd6881db36b5e08329cf23cecdd vs subject 78975d5144a1cd12e98898d573cf6536
    #
    # Expected matrices from pyani v0.2 give us expected values of:
    # identity 0.9945938461538463, aln_length 39169, sim_errors 215
    # Note final digit wobble from 3 (old) to 2 (new)


def test_running_anib(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check can compute and log column of ANIb comparisons to DB."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "new.sqlite"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "anib.json"

    tool = tools.get_blastn()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus anib ...",
        status="Testing",
        name="Testing logging_anib",
        method="ANIb",
        program=tool.exe_path.stem,
        version=tool.version,
        fragsize=anib.FRAGSIZE,  # will be used by default in log_anib
        create_db=True,
    )

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = session.query(db_orm.Run).one()
    assert run.run_id == 1
    hash_to_filename = {_.genome_hash: _.fasta_filename for _ in run.fasta_hashes}
    hash_to_length = {_.genome_hash: _.length for _ in run.genomes}

    subject_hash = list(hash_to_filename)[1]
    private_cli.compute_anib(
        logger,
        tmp_dir,
        session,
        run,
        tmp_json,
        input_genomes_tiny,
        hash_to_filename,
        {},  # not used for ANIb
        query_hashes=hash_to_length,  # order should not matter!
        subject_hash=subject_hash,
    )

    private_cli.import_json_comparisons(logger, session, tmp_json)

    assert session.query(db_orm.Comparison).count() == 3  # noqa: PLR2004
    assert (
        session.query(db_orm.Comparison)
        .where(db_orm.Comparison.subject_hash == subject_hash)
        .count()
        == 3  # noqa: PLR2004
    )

    # Check the intermediate fragmented FASTA files match
    for fname in (input_genomes_tiny / "intermediates/ANIb").glob("*-fragments.fna"):
        # Intermediate uses f"{query_hash}-fragments-{fragsize}-pid{os.getpid()}.fna"
        tmp_frag_files = list(
            tmp_dir.glob(f"{fname.name.rsplit('-', 1)[0]}-fragments-*.fna")
        )
        # Only computed one row, so should be only one copy of the fragment file
        assert len(tmp_frag_files) == 1, (tmp_frag_files, fname.name)
        assert filecmp.cmp(fname, tmp_frag_files[0])

    # Check the intermediate TSV files from blastn match
    for fname in (input_genomes_tiny / "intermediates/ANIb").glob(
        f"*_vs_{subject_hash}.tsv"
    ):
        assert filecmp.cmp(fname, tmp_dir / fname.name)

    # No real need to test the ANI values here, will be done elsewhere.
    for query_hash, query_filename in hash_to_filename.items():
        pytest.approx(
            get_matrix_entry(
                input_genomes_tiny / "matrices/ANIb_identity.tsv",
                Path(query_filename).stem,
                Path(hash_to_filename[subject_hash]).stem,
            )
            == session.query(db_orm.Comparison)
            .where(db_orm.Comparison.query_hash == query_hash)
            .where(db_orm.Comparison.subject_hash == subject_hash)
            .one()
            .identity
        )
    session.close()
    tmp_db.unlink()
