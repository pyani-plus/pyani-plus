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
"""Test methods for calculating ANI using minimap2.

These tests are intended to be run from the repository root using:

make test
"""

from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, setup_logger, tools
from pyani_plus.methods import animinimap2


def test_animinimap2_parsing_empty(tmp_path: str) -> None:
    """Check parsing of empty animinimap2 output files."""
    tmp_empty = Path(tmp_path) / "empty.tsv"
    with tmp_empty.open("wb"):
        pass

    with pytest.raises(ValueError, match="Empty PAF file from minimap2"):
        animinimap2.parse_minimap2_paf_file(tmp_empty)


def test_running_animinimap2(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check that animinimap2 can be run on test input genomes."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "new.sqlite"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "animinimap2.json"

    tool = tools.get_minimap2()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus animinimap2 ...",
        status="Testing",
        name="Testing animinimap2",
        method="animinimap2",
        program=tool.exe_path.stem,
        version=tool.version,
        create_db=True,
        mode=animinimap2.DEFAULT_PRESET,
    )

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = session.query(db_orm.Run).one()
    assert run.run_id == 1
    filename_to_hash = {_.fasta_filename: _.genome_hash for _ in run.fasta_hashes}
    hash_to_filename = {_.genome_hash: _.fasta_filename for _ in run.fasta_hashes}
    hash_to_lengths = {_.genome_hash: _.length for _ in run.genomes}

    private_cli.compute_animinimap2(
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
