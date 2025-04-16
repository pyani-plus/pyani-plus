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
"""Tests for the sourmash implementation.

These tests are intended to be run from the repository root using:

make test
"""

# Required to support pytest automated testing
from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, setup_logger, tools
from pyani_plus.methods import sourmash


def test_prepare_genomes_bad_method(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check error handling of sourmash.prepare_genomes with wrong method."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad-args.db"

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus sourmash ...",
        status="Testing",
        name="Testing sourmash prepare-genomes",
        method="guessing",
        program="guestimator",
        version="0.0a1",
        create_db=True,
    )
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.load_run(session, run_id=1)
    with pytest.raises(
        SystemExit,
        match="Expected run to be for sourmash, not method guessing",
    ):
        next(
            sourmash.prepare_genomes(logger, run, tmp_dir)
        )  # should error before checks cache
    session.close()


def test_prepare_genomes_bad_kmer(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check error handling of sourmash.prepare_genomes without k-mer size."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad-args.db"

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus sourmash ...",
        status="Testing",
        name="Testing sourmash prepare-genomes",
        method="sourmash",
        program="sourmash",
        version="0.0a1",
        extra="scaled=" + str(sourmash.SCALED),
        create_db=True,
    )
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.load_run(session, run_id=1)
    with pytest.raises(
        SystemExit,
        match=f"sourmash requires a k-mer size, default is {sourmash.KMER_SIZE}",
    ):
        next(
            sourmash.prepare_genomes(logger, run, cache=tmp_dir)
        )  # should error before checks cache
    session.close()


def test_prepare_genomes_bad_cache(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check error handling of sourmash.prepare_genomes without k-mer size."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad-args.db"

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus sourmash ...",
        status="Testing",
        name="Testing sourmash prepare-genomes",
        method="sourmash",
        program="sourmash",
        version="0.0a1",
        kmersize=sourmash.KMER_SIZE,
        extra="scaled=" + str(sourmash.SCALED),
        create_db=True,
    )
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.load_run(session, run_id=1)
    with pytest.raises(
        ValueError,
        match="Cache directory '/does/not/exist' does not exist",
    ):
        next(
            sourmash.prepare_genomes(logger, run, cache=Path("/does/not/exist"))
        )  # should error before checks cache
    session.close()


def test_prepare_genomes_bad_extra(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check error handling of sourmash.prepare_genomes without extra (scaling)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad-args.db"

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus sourmash ...",
        status="Testing",
        name="Testing sourmash prepare-genomes",
        method="sourmash",
        program="sourmash",
        version="0.0a1",
        kmersize=sourmash.KMER_SIZE,
        create_db=True,
    )
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.load_run(session, run_id=1)
    with pytest.raises(
        SystemExit,
        match=f"sourmash requires extra setting, default is scaled={sourmash.SCALED}",
    ):
        next(
            sourmash.prepare_genomes(logger, run, cache=tmp_dir)
        )  # should error before checks cache
    session.close()


def test_parser_with_bad_branchwater(tmp_path: str) -> None:
    """Check self-vs-self is one in sourmash compare parser."""
    mock_csv = Path(tmp_path) / "faked.csv"
    with mock_csv.open("w") as handle:
        handle.write(
            "max_containment_ani,query_name,match_name,query_containment_ani\n"
        )
        handle.write("\n")  # parser will skip blank lines
        handle.write("1.0,AAAAAA,AAAAAA,1.0\n")
        handle.write("0.9,AAAAAA,BBBBBB,0.85\n")
        handle.write("NaN,BBBBBB,BBBBBB,NaN\n")  # fails self-vs-self 100%
    expected = {
        ("AAAAAA", "AAAAAA"),
        ("AAAAAA", "BBBBBB"),
        ("BBBBBB", "AAAAAA"),
        ("BBBBBB", "BBBBBB"),
    }
    logger = setup_logger(None)
    parser = sourmash.parse_sourmash_manysearch_csv(logger, mock_csv, expected)
    assert next(parser) == ("AAAAAA", "AAAAAA", 1.0, 1.0)
    assert next(parser) == ("AAAAAA", "BBBBBB", 0.85, 0.9)
    with pytest.raises(
        ValueError,
        match="Expected sourmash manysearch BBBBBB vs self to be one, not 'NaN'",
    ):
        next(parser)

    # Now tell it just expect one entry...
    parser = sourmash.parse_sourmash_manysearch_csv(
        logger, mock_csv, {("AAAAAA", "AAAAAA")}
    )
    assert next(parser) == ("AAAAAA", "AAAAAA", 1.0, 1.0)
    with pytest.raises(
        SystemExit, match="Did not expect AAAAAA vs BBBBBB in faked.csv"
    ):
        next(parser)


def test_parser_with_bad_header(tmp_path: str) -> None:
    """Check sourmash branchwater parser with bad header."""
    mock_csv = Path(tmp_path) / "faked.csv"
    with mock_csv.open("w") as handle:
        # Sourmash branchwater does not use subject_containment_ani,
        # rather they have query_containment_ani and match_containment_ani
        handle.write(
            "max_containment_ani,query_name,match_name,subject_containment_ani\n"
        )
    logger = setup_logger(None)
    parser = sourmash.parse_sourmash_manysearch_csv(logger, mock_csv, set())
    with pytest.raises(
        SystemExit,
        match="Missing expected fields in sourmash manysearch header, found: "
        "'max_containment_ani,query_name,match_name,subject_containment_ani'",
    ):
        next(parser)


def test_compute_bad_args(tmp_path: str) -> None:
    """Check compute_sourmash error handling."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad args.db"
    tmp_json = tmp_dir / "bad args.json"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.Run()  # empty
    tool = tools.get_sourmash()
    config = db_orm.Configuration(
        method="sourmash",
        program=tool.exe_path.name,
        version=tool.version,
        kmersize=31,
        extra="scaled=1234",
    )
    run = db_orm.Run(configuration=config)
    with pytest.raises(
        SystemExit,
        match=(
            "Missing sourmash signatures directory"
            f" '{tmp_dir}/sourmash_k=31_scaled=1234' - check cache setting."
        ),
    ):
        private_cli.compute_sourmash(
            logger,
            tmp_dir,
            session,
            run,
            tmp_json,
            tmp_dir,
            {},
            {},
            {"ABCDE": 12345},
            "HIJKL",
            cache=tmp_dir,
        )


def test_compute_tile_bad_args(tmp_path: str) -> None:
    """Check compute_sourmash_tile error handling."""
    tmp_dir = Path(tmp_path)
    tool = tools.ExternalToolData(exe_path=Path("sourmash"), version="0.0a1")
    logger = setup_logger(None)
    with pytest.raises(
        ValueError, match="Given cache directory '/does/not/exist' does not exist"
    ):
        next(
            sourmash.compute_sourmash_tile(
                logger,
                tool,
                {
                    "",
                },
                {
                    "",
                },
                Path("/does/not/exist"),
                tmp_dir,
            )
        )
    with pytest.raises(SystemExit, match="Return code 1 from: sourmash sig collect "):
        next(
            sourmash.compute_sourmash_tile(
                logger,
                tool,
                {
                    "ACBDE",
                },
                {
                    "ABCDE",
                },
                tmp_dir,
                tmp_dir,
            )
        )


def test_compute_tile_stale_cvs(
    caplog: pytest.LogCaptureFixture, tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check compute_sourmash_tile with stale sig-lists."""
    tmp_dir = Path(tmp_path)

    query_csv = tmp_dir / "query_sigs.csv"
    query_csv.touch()
    subject_csv = tmp_dir / "subject_sigs.csv"
    subject_csv.touch()

    tool = tools.get_sourmash()
    logger = setup_logger(None)
    next(
        sourmash.compute_sourmash_tile(
            logger,
            tool,
            {"689d3fd6881db36b5e08329cf23cecdd", "5584c7029328dc48d33f95f0a78f7e57"},
            {"689d3fd6881db36b5e08329cf23cecdd", "78975d5144a1cd12e98898d573cf6536"},
            input_genomes_tiny / "intermediates/sourmash",
            tmp_dir,
        )
    )
    output = caplog.text
    assert f"Race condition? Replacing intermediate file '{query_csv}'" in output, (
        output
    )
    assert f"Race condition? Replacing intermediate file '{subject_csv}'" in output, (
        output
    )
