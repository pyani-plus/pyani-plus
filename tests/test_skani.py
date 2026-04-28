# The MIT License
#
# Copyright (c) 2024-2026 University of Strathclyde
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
"""Test methods for calculating ANI using skani.

These tests are intended to be run from the repository root using:

make test
"""

from pathlib import Path

from pyani_plus import db_orm, private_cli, setup_logger, tools
from pyani_plus.methods import skani


def test_skani_parsing(input_genomes_tiny: Path) -> None:
    """Check parsing of test skani output files."""
    assert skani.parse_skani(
        input_genomes_tiny
        / "intermediates"
        / "skani"
        / "689d3fd6881db36b5e08329cf23cecdd_vs_689d3fd6881db36b5e08329cf23cecdd.skani"
    ) == (
        "../fixtures/viral_example/MGV-GENOME-0264574.fas",
        "../fixtures/viral_example/MGV-GENOME-0264574.fas",
        0.01 * 100.00,
        0.01 * 99.65,
        0.01 * 99.65,
        "MGV_MGV-GENOME-0264574",
        "MGV_MGV-GENOME-0264574",
    )

    assert skani.parse_skani(
        input_genomes_tiny
        / "intermediates"
        / "skani"
        / "78975d5144a1cd12e98898d573cf6536_vs_5584c7029328dc48d33f95f0a78f7e57.skani"
    ) == (
        "../fixtures/viral_example/OP073605.fasta",
        "../fixtures/viral_example/MGV-GENOME-0266457.fna",
        0.01 * 99.64,
        0.01 * 68.28,
        0.01 * 99.66,
        "OP073605.1 MAG: Bacteriophage sp. isolate 0984_12761, complete genome",
        "MGV_MGV-GENOME-0266457",
    )


def test_skani_parsing_empty(input_gzip_bacteria: Path) -> None:
    """Check parsing of empty skani output files."""
    assert skani.parse_skani(
        input_gzip_bacteria
        / "intermediates"
        / "skani"
        / "9a9e23bfc5a184b8149e07e267d133b0_vs_073194224aa8c13bebc1d14a3e74a3e7.skani"
    ) == ("", "", None, None, None, "", "")


def test_running_skani(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check that skani can be run on test input genomes."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "new.sqlite"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "skani.json"

    tool = tools.get_skani()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus skani ...",
        status="Testing",
        name="Testing skani",
        method="skani",
        program=tool.exe_path.stem,
        version=tool.version,
        create_db=True,
        mode=skani.MODE,
    )

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = session.query(db_orm.Run).one()
    assert run.run_id == 1
    filename_to_hash = {_.fasta_filename: _.genome_hash for _ in run.fasta_hashes}
    hash_to_filename = {_.genome_hash: _.fasta_filename for _ in run.fasta_hashes}
    hash_to_lengths = {_.genome_hash: _.length for _ in run.genomes}

    private_cli.compute_skani(
        logger,
        tmp_dir,
        session,
        run,
        tmp_json,
        input_genomes_tiny,
        hash_to_filename,
        filename_to_hash,
        query_hashes=hash_to_lengths,
        subject_hash=list(hash_to_filename)[1],
    )


def test_running_skani_gzip(
    tmp_path: str,
    input_gzip_bacteria: Path,
) -> None:
    """Check that skani can be run on gzip test input genomes."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "new.sqlite"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "skani.json"

    tool = tools.get_skani()

    private_cli.log_run(
        fasta=input_gzip_bacteria,
        database=tmp_db,
        cmdline="pyani-plus skani ...",
        status="Testing",
        name="Testing skani",
        method="skani",
        program=tool.exe_path.stem,
        version=tool.version,
        create_db=True,
        mode=skani.MODE,
    )

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = session.query(db_orm.Run).one()
    assert run.run_id == 1
    filename_to_hash = {_.fasta_filename: _.genome_hash for _ in run.fasta_hashes}
    hash_to_filename = {_.genome_hash: _.fasta_filename for _ in run.fasta_hashes}
    hash_to_lengths = {_.genome_hash: _.length for _ in run.genomes}

    private_cli.compute_skani(
        logger,
        tmp_dir,
        session,
        run,
        tmp_json,
        input_gzip_bacteria,
        hash_to_filename,
        filename_to_hash,
        query_hashes=hash_to_lengths,
        subject_hash=list(hash_to_filename)[1],
    )
