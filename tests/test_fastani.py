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
"""Tests for the fastANI implementation.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, tools
from pyani_plus.methods.method_fastani import parse_fastani_file


def test_bad_path() -> None:
    """Confirm giving an empty path etc fails."""
    with pytest.raises(FileNotFoundError, match="No such file or directory:"):
        list(parse_fastani_file(Path("/does/not/exist"), {}))


def test_empty_path() -> None:
    """Confirm giving an empty path etc fails."""
    with pytest.raises(ValueError, match="Input file /dev/null is empty"):
        list(parse_fastani_file(Path("/dev/null"), {}))


def test_missing_db(tmp_path: str) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.fastani(
            database=tmp_db,
            run_id=1,
            subject="XXX",
        )


def test_running_fastani(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check can run a fastANI comparison to DB."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_fastani()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus fastani ...",
        status="Testing",
        name="Testing log_fastani",
        method="fastANI",
        program=tool.exe_path.stem,
        version=tool.version,
        fragsize=1000,
        kmersize=13,  # must be at most 16
        minmatch=0.9,
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    with pytest.raises(
        SystemExit,
        match="ERROR: Did not recognise 'XXXX' as an MD5 hash or filename in run-id 1",
    ):
        private_cli.fastani(
            database=tmp_db,
            run_id=1,
            subject="XXXX",
        )

    private_cli.fastani(
        database=tmp_db,
        run_id=1,
        subject="MGV-GENOME-0266457.fna",  # will test using a hash next
    )

    # Check the recorded comparison values
    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 3  # noqa: PLR2004
    # No need to test the ANI values here, will be done elsewhere.

    # Do another row, should accept a hash:
    private_cli.fastani(
        database=tmp_db, run_id=1, subject="689d3fd6881db36b5e08329cf23cecdd"
    )
    assert session.query(db_orm.Comparison).count() == 6  # noqa: PLR2004

    # Do the same row again, should skip gracefully:
    private_cli.fastani(
        database=tmp_db, run_id=1, subject="689d3fd6881db36b5e08329cf23cecdd"
    )
    assert session.query(db_orm.Comparison).count() == 6  # noqa: PLR2004

    session.close()
    tmp_db.unlink()
