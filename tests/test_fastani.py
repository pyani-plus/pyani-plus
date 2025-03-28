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
"""Tests for the fastANI implementation.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, setup_logger, tools
from pyani_plus.methods.fastani import parse_fastani_file


def test_bad_path() -> None:
    """Confirm giving an empty path etc fails."""
    with pytest.raises(FileNotFoundError, match="No such file or directory:"):
        list(parse_fastani_file(Path("/does/not/exist"), {}, set()))


def test_running_fastani(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check can run a fastANI comparison to DB."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "new.sqlite"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "fastani.json"

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

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = session.query(db_orm.Run).one()
    assert run.run_id == 1
    filename_to_hash = {_.fasta_filename: _.genome_hash for _ in run.fasta_hashes}
    hash_to_filename = {_.genome_hash: _.fasta_filename for _ in run.fasta_hashes}
    hash_to_lengths = {_.genome_hash: _.length for _ in run.genomes}

    private_cli.compute_fastani(
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

    private_cli.import_json_comparisons(logger, session, tmp_json)

    assert session.query(db_orm.Comparison).count() == 3  # noqa: PLR2004

    session.close()
    tmp_db.unlink()
