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

from pyani_plus import db_orm, private_cli
from pyani_plus.methods import sourmash
from pyani_plus.tools import get_sourmash


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


def test_missing_db(tmp_path: str) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.compute_tile(
            database=tmp_db,
            run_id=1,
            tile=0,
        )


def test_sketch_signatures_bad_args(tmp_path: str) -> None:
    """Try calling sketch_signatures with bad configuration."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad.db"

    tool = get_sourmash()

    private_cli.log_configuration(
        database=tmp_db,
        method="sourmash",
        program=tool.exe_path.stem,
        version=tool.version,
        create_db=True,
    )

    session = db_orm.connect_to_db(tmp_db)
    config = session.query(db_orm.Configuration).one()
    assert not config.kmersize
    assert not config.extra
    with pytest.raises(
        SystemExit, match="ERROR: sourmash requires a k-mer size, default is 31"
    ):
        next(sourmash.sketch_signatures(tool, config, {}, tmp_dir))

    private_cli.log_configuration(
        database=tmp_db,
        method="sourmash",
        program=tool.exe_path.stem,
        version=tool.version,
        kmersize=17,
    )

    session = db_orm.connect_to_db(tmp_db)
    config = (
        session.query(db_orm.Configuration)
        .where(db_orm.Configuration.configuration_id == 2)  # noqa: PLR2004
        .one()
    )
    assert config.kmersize == 17  # noqa: PLR2004
    assert not config.extra
    with pytest.raises(
        SystemExit,
        match="ERROR: sourmash requires scaled or num, default is scaled=1000",
    ):
        next(sourmash.sketch_signatures(tool, config, {}, tmp_dir))


def test_compute_with_wrong_version(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check mismatched sourmash version fails."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus sourmash ...",
        status="Testing",
        name="Testing log_sourmash",
        method="sourmash",
        program="sourmash",
        version="42",
        kmersize=sourmash.KMER_SIZE,
        extra="scaled=" + str(sourmash.SCALED),
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    with pytest.raises(
        SystemExit,
        match="ERROR: Run configuration was sourmash 42 but we have sourmash 4.",
    ):
        private_cli.compute_tile(
            database=tmp_db,
            run_id=1,
            tile=0,
        )
