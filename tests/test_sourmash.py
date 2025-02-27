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

from pyani_plus import db_orm, private_cli, tools
from pyani_plus.methods import sourmash


def test_prepare_genomes_bad_method(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
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
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")
    session = db_orm.connect_to_db(tmp_db)
    run = db_orm.load_run(session, run_id=1)

    with pytest.raises(
        SystemExit,
        match="ERROR: Expected run to be for sourmash, not method guessing",
    ):
        next(sourmash.prepare_genomes(run, tmp_dir))  # should error before checks cache
    session.close()


def test_prepare_genomes_bad_kmer(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
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
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")
    session = db_orm.connect_to_db(tmp_db)
    run = db_orm.load_run(session, run_id=1)

    with pytest.raises(
        SystemExit,
        match=f"ERROR: sourmash requires a k-mer size, default is {sourmash.KMER_SIZE}",
    ):
        next(
            sourmash.prepare_genomes(run, cache=tmp_dir)
        )  # should error before checks cache
    session.close()


def test_prepare_genomes_bad_cache(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
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
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")
    session = db_orm.connect_to_db(tmp_db)
    run = db_orm.load_run(session, run_id=1)

    with pytest.raises(
        ValueError,
        match="ERROR: Cache directory /does/not/exist does not exist",
    ):
        next(
            sourmash.prepare_genomes(run, cache=Path("/does/not/exist"))
        )  # should error before checks cache
    session.close()


def test_prepare_genomes_bad_extra(
    capsys: pytest.CaptureFixture[str],
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
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")
    session = db_orm.connect_to_db(tmp_db)
    run = db_orm.load_run(session, run_id=1)

    with pytest.raises(
        SystemExit,
        match=f"ERROR: sourmash requires extra setting, default is scaled={sourmash.SCALED}",
    ):
        next(
            sourmash.prepare_genomes(run, cache=tmp_dir)
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
    parser = sourmash.parse_sourmash_manysearch_csv(mock_csv, expected)
    assert next(parser) == ("AAAAAA", "AAAAAA", 1.0, 1.0)
    assert next(parser) == ("AAAAAA", "BBBBBB", 0.85, 0.9)
    with pytest.raises(
        ValueError,
        match="Expected sourmash manysearch BBBBBB vs self to be one, not 'NaN'",
    ):
        next(parser)

    # Now tell it just expect one entry...
    parser = sourmash.parse_sourmash_manysearch_csv(mock_csv, {("AAAAAA", "AAAAAA")})
    assert next(parser) == ("AAAAAA", "AAAAAA", 1.0, 1.0)
    with pytest.raises(
        ValueError, match="Did not expect AAAAAA vs BBBBBB in faked.csv"
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
    parser = sourmash.parse_sourmash_manysearch_csv(mock_csv, set())
    with pytest.raises(
        SystemExit,
        match="ERROR - Missing expected fields in sourmash manysearch header, found: "
        "'max_containment_ani,query_name,match_name,subject_containment_ani'",
    ):
        next(parser)


def test_compute_bad_args(tmp_path: str) -> None:
    """Check compute_sourmash error handling."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad_args.db"
    tool = tools.ExternalToolData(exe_path=Path("sourmash"), version="0.0a1")
    session = db_orm.connect_to_db(tmp_db)
    run = db_orm.Run()  # empty
    with pytest.raises(SystemExit, match="ERROR: Not given a cache directory"):
        private_cli.compute_sourmash(
            tmp_dir, session, run, tmp_dir, {}, {}, {"ABCDE": 12345}, "HIJKL"
        )
    with pytest.raises(
        SystemExit,
        match="ERROR: Cache directory /does/not/exist does not exist - check cache setting.",
    ):
        private_cli.compute_sourmash(
            tmp_dir,
            session,
            run,
            tmp_dir,
            {},
            {},
            {"ABCDE": 12345},
            "HIJKL",
            cache=Path("/does/not/exist"),
        )

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
            "ERROR: Missing sourmash signatures directory"
            f" {tmp_dir}/sourmash_k=31_scaled=1234 - check cache setting."
        ),
    ):
        private_cli.compute_sourmash(
            tmp_dir,
            session,
            run,
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
    with pytest.raises(
        ValueError, match="Given cache directory /does/not/exist does not exist"
    ):
        next(
            sourmash.compute_sourmash_tile(
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
    with pytest.raises(
        SystemExit, match="ERROR: Return code 1 from: sourmash sig collect "
    ):
        next(
            sourmash.compute_sourmash_tile(
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
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check compute_sourmash_tile with stale sig-lists."""
    tmp_dir = Path(tmp_path)

    query_csv = tmp_dir / "query_sigs.csv"
    query_csv.touch()
    subject_csv = tmp_dir / "subject_sigs.csv"
    subject_csv.touch()

    tool = tools.get_sourmash()
    next(
        sourmash.compute_sourmash_tile(
            tool,
            {"689d3fd6881db36b5e08329cf23cecdd", "5584c7029328dc48d33f95f0a78f7e57"},
            {"689d3fd6881db36b5e08329cf23cecdd", "78975d5144a1cd12e98898d573cf6536"},
            input_genomes_tiny / "intermediates/sourmash",
            tmp_dir,
        )
    )
    output = capsys.readouterr().err
    assert (
        f"WARNING: Race condition? Replacing intermediate file {query_csv}" in output
    ), output
    assert (
        f"WARNING: Race condition? Replacing intermediate file {subject_csv}" in output
    ), output
