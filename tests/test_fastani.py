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

from . import get_matrix_entry


def test_bad_path() -> None:
    """Confirm giving an empty path etc fails."""
    with pytest.raises(FileNotFoundError, match="No such file or directory:"):
        list(parse_fastani_file(Path("/does/not/exist"), {}))


def test_empty_path() -> None:
    """Confirm giving an empty path etc fails."""
    with pytest.raises(ValueError, match="Input file /dev/null is empty"):
        list(parse_fastani_file(Path("/dev/null"), {}))


def test_missing_db(tmp_path: str, fastani_targets_indir: Path) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.log_fastani(
            database=tmp_db,
            run_id=1,
            fastani=fastani_targets_indir / "all_vs_MGV-GENOME-0266457.fastani",
        )


def test_logging_fastani(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
    fastani_targets_indir: Path,
) -> None:
    """Check can log a fastANI comparison to DB."""
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
        kmersize=51,
        minmatch=0.9,
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    private_cli.log_fastani(
        database=tmp_db,
        run_id=1,
        fastani=fastani_targets_indir / "all_vs_MGV-GENOME-0266457.fastani",
    )

    # Check the recorded comparison values
    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 3  # noqa: PLR2004
    query = "689d3fd6881db36b5e08329cf23cecdd"  # MGV-GENOME-0264574.fas
    subject = "78975d5144a1cd12e98898d573cf6536"  # MGV-GENOME-0266457.fna
    comp = (
        session.query(db_orm.Comparison)
        .where(db_orm.Comparison.query_hash == query)
        .one()
    )
    assert comp.subject_hash == subject
    assert comp.configuration_id == 1  # by construction

    # Using approx to avoid 0.9950140000000001 != 0.995014
    pytest.approx(
        comp.identity,
        get_matrix_entry(
            input_genomes_tiny / "matrices" / "fastANI_identity.tsv", query, subject
        ),
    )
    pytest.approx(
        comp.aln_length,
        get_matrix_entry(
            input_genomes_tiny / "matrices" / "fastANI_aln_lengths.tsv", query, subject
        ),
    )
    pytest.approx(
        comp.sim_errors,
        get_matrix_entry(
            input_genomes_tiny / "matrices" / "fastANI_sim_errors.tsv", query, subject
        ),
    )
    pytest.approx(
        comp.cov_query,
        get_matrix_entry(
            input_genomes_tiny / "matrices" / "fastANI_coverage.tsv", query, subject
        ),
    )
    session.close()
    tmp_db.unlink()
