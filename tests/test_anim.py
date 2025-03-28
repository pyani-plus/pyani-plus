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
"""Test methods for calculating ANIm.

These tests are intended to be run from the repository root using:

make test
"""

import filecmp
from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, setup_logger, tools, utils
from pyani_plus.methods import anim

from . import get_matrix_entry


@pytest.fixture
def aligned_regions() -> dict:
    """Example of aligned regions with overlaps."""  # noqa: D401
    return {
        "MGV-GENOME-0266457": [(1, 37636), (17626, 39176)],
    }


@pytest.fixture
def genome_hashes(input_genomes_tiny: Path) -> dict:
    """Return the MD5 checksum hash digest of the input genomes."""
    return {
        record.stem: utils.file_md5sum(record)
        for record in input_genomes_tiny.iterdir()
    }


def test_delta_parsing(input_genomes_tiny: Path) -> None:
    """Check parsing of test NUCmer .delta/.filter file."""
    assert anim.parse_delta(
        input_genomes_tiny
        / "intermediates/ANIm/689d3fd6881db36b5e08329cf23cecdd_vs_78975d5144a1cd12e98898d573cf6536.filter",
    ) == (
        39169,
        39176,
        0.9962487643734,
        222,
    )

    with pytest.raises(ValueError, match="Empty delta file from nucmer, /dev/null"):
        anim.parse_delta(Path("/dev/null"))


def test_bad_alignments_parsing(input_genomes_bad_alignments: Path) -> None:
    """Check parsing of test NUCmer .delta/.filter file."""
    assert anim.parse_delta(
        input_genomes_bad_alignments
        / "intermediates/ANIm/689d3fd6881db36b5e08329cf23cecdd_vs_a30481565b45f6bbc6ce5260503067e0.filter",
    ) == (None, None, None, None)


def test_aligned_bases_count(aligned_regions: dict) -> None:
    """Check only aligned bases in non-overlapping regions are counted."""
    assert anim.get_aligned_bases_count(aligned_regions) == 39176  # noqa: PLR2004


def test_running_anim(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check can compute and log column of ANIm comparisons to DB."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "new.sqlite"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "anim.json"

    tool = tools.get_nucmer()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus anim ...",
        status="Testing",
        name="Testing anim",
        method="ANIm",
        program=tool.exe_path.stem,
        version=tool.version,
        mode=anim.MODE,
        create_db=True,
    )

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = session.query(db_orm.Run).one()
    assert run.run_id == 1
    hash_to_filename = {_.genome_hash: _.fasta_filename for _ in run.fasta_hashes}
    hash_to_length = {_.genome_hash: _.length for _ in run.genomes}

    subject_hash = list(hash_to_filename)[1]
    private_cli.compute_anim(
        logger,
        tmp_dir,
        session,
        run,
        tmp_json,
        input_genomes_tiny,
        hash_to_filename,
        {},  # not used for ANIm
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

    # Could check nucmer output against target fixtures?

    # Check the intermediate delta-filter files
    subject_stem = Path(hash_to_filename[subject_hash]).stem
    for fname in (input_genomes_tiny / "intermediates/ANIb").glob(
        f"*_vs_{subject_stem}.delta"
    ):
        assert (tmp_dir / fname.name).is_file(), list(tmp_dir.glob("*"))
        assert filecmp.cmp(fname, tmp_dir / fname.name)

    # No real need to test the ANI values here, will be done elsewhere.
    for query_hash, query_filename in hash_to_filename.items():
        pytest.approx(
            get_matrix_entry(
                input_genomes_tiny / "matrices/ANIm_identity.tsv",
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
