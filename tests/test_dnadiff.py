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
"""Test methods for calculating dnadiff.

These tests are intended to be run from the repository root using:

make test
"""

# Required to support pytest automated testing
from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, setup_logger, tools
from pyani_plus.methods import dnadiff

from . import get_matrix_entry


@pytest.fixture
def expected_mcoords_output() -> tuple[float, int]:
    """Expected average identity and aligned bases with gaps.

    MGV-GENOME-0264574 (REF) MGV-GENOME-0266457 (QRY).
    """  # noqa: D401
    return (0.996266174669021, 39253)


@pytest.fixture
def expected_gap_lengths_qry() -> int:
    """Expected gap length in the QRY alignment."""  # noqa: D401
    return 84


def test_parse_mcoords(
    input_genomes_tiny: Path, expected_mcoords_output: tuple[float, int]
) -> None:
    """Check parsing of test mcoords file."""
    assert expected_mcoords_output == dnadiff.parse_mcoords(
        input_genomes_tiny
        / "intermediates/dnadiff/689d3fd6881db36b5e08329cf23cecdd_vs_78975d5144a1cd12e98898d573cf6536.mcoords"
    )


def test_parse_mcoords_bad_alignment(input_genomes_bad_alignments: Path) -> None:
    """Check parsing of test mcoords file for bad alignments example."""
    assert dnadiff.parse_mcoords(
        input_genomes_bad_alignments
        / "intermediates/dnadiff/689d3fd6881db36b5e08329cf23cecdd_vs_a30481565b45f6bbc6ce5260503067e0.mcoords"
    ) == (None, None)


def test_parse_qdiff(input_genomes_tiny: Path, expected_gap_lengths_qry: int) -> None:
    """Check parsing of test qdiff file."""
    assert expected_gap_lengths_qry == dnadiff.parse_qdiff(
        input_genomes_tiny
        / "intermediates/dnadiff/689d3fd6881db36b5e08329cf23cecdd_vs_78975d5144a1cd12e98898d573cf6536.qdiff"
    )


def test_parse_qdiff_bad_alignments(input_genomes_bad_alignments: Path) -> None:
    """Check parsing of test qdiff file for bad alignments example."""
    assert (
        dnadiff.parse_qdiff(
            input_genomes_bad_alignments
            / "intermediates/dnadiff/689d3fd6881db36b5e08329cf23cecdd_vs_a30481565b45f6bbc6ce5260503067e0.qdiff"
        )
        is None
    )


def test_running_dnadiff(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check can compute and log column of dnadiff comparisons to DB."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "new.sqlite"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "dna differences.json"

    tool = tools.get_nucmer()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus dnadiff ...",
        status="Testing",
        name="Testing dnadiff",
        method="dnadiff",
        program=tool.exe_path.stem,
        version=tool.version,
        create_db=True,
    )

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = session.query(db_orm.Run).one()
    assert run.run_id == 1
    hash_to_filename = {_.genome_hash: _.fasta_filename for _ in run.fasta_hashes}
    hash_to_length = {_.genome_hash: _.length for _ in run.genomes}

    subject_hash = list(hash_to_filename)[1]
    private_cli.compute_dnadiff(
        logger,
        tmp_dir,
        session,
        run,
        tmp_json,
        input_genomes_tiny,
        hash_to_filename,
        {},  # not used for dnadiff
        query_hashes=hash_to_length,  # order should not matter!
        subject_hash=subject_hash,
    )
    assert tmp_json.is_file()

    private_cli.import_json_comparisons(logger, session, tmp_json)

    assert session.query(db_orm.Comparison).count() == 3  # noqa: PLR2004
    assert (
        session.query(db_orm.Comparison)
        .where(db_orm.Comparison.subject_hash == subject_hash)
        .count()
        == 3  # noqa: PLR2004
    )

    # Check the intermediate files match?

    # No real need to test the ANI values here, will be done elsewhere.
    for query_hash, query_filename in hash_to_filename.items():
        pytest.approx(
            get_matrix_entry(
                input_genomes_tiny / "matrices/dnadiff_identity.tsv",
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
