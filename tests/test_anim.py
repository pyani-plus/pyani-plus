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
"""Test methods for calculating ANIm.

These tests are intended to be run from the repository root using:

make test
"""

# Required to support pytest automated testing
from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, tools, utils
from pyani_plus.methods import method_anim

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
    assert method_anim.parse_delta(
        input_genomes_tiny
        / "intermediates/ANIm/MGV-GENOME-0264574_vs_MGV-GENOME-0266457.filter",
    ) == (
        39169,
        39176,
        0.9962487643734,
        222,
    )

    with pytest.raises(ValueError, match="Empty delta file from nucmer, /dev/null"):
        method_anim.parse_delta(Path("/dev/null"))


def test_bad_alignments_parsing(input_genomes_bad_alignments: Path) -> None:
    """Check parsing of test NUCmer .delta/.filter file."""
    assert method_anim.parse_delta(
        input_genomes_bad_alignments
        / "intermediates/ANIm/MGV-GENOME-0264574_vs_MGV-GENOME-0357962.filter",
    ) == (None, None, None, None)


def test_aligned_bases_count(aligned_regions: dict) -> None:
    """Check only aligned bases in non-overlapping regions are counted."""
    assert method_anim.get_aligned_bases_count(aligned_regions) == 39176  # noqa: PLR2004


def test_missing_db(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.log_anim(
            database=tmp_db,
            run_id=1,
            # These are for the comparison table
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            deltafilter=input_genomes_tiny
            / "intermediates/ANIm/MGV-GENOME-0264574_vs_MGV-GENOME-0266457.filter",
        )


def test_bad_query_or_subject(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Mismatch between query or subject FASTA in fastANI output and commandline."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --query-fasta .*/MGV-GENOME-0266457.fna"
            " but query in deltafilter filename was MGV-GENOME-0264574"
        ),
    ):
        private_cli.log_anim(
            database=tmp_db,
            run_id=1,
            query_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            deltafilter=input_genomes_tiny
            / "intermediates/ANIm/MGV-GENOME-0264574_vs_MGV-GENOME-0266457.filter",
        )

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --subject-fasta .*/MGV-GENOME-0264574.fas"
            " but subject in deltafilter filename was MGV-GENOME-0266457"
        ),
    ):
        private_cli.log_anim(
            database=tmp_db,
            run_id=1,
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            deltafilter=input_genomes_tiny
            / "intermediates/ANIm/MGV-GENOME-0264574_vs_MGV-GENOME-0266457.filter",
        )


def test_logging_anim(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check can log a ANIm comparison to DB."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_nucmer()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus anim ...",
        status="Testing",
        name="Testing log_anim",
        method="ANIm",
        program=tool.exe_path.stem,
        version=tool.version,
        mode=method_anim.MODE,
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    private_cli.log_anim(
        database=tmp_db,
        run_id=1,
        query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
        subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
        deltafilter=input_genomes_tiny
        / "intermediates/ANIm/MGV-GENOME-0264574_vs_MGV-GENOME-0266457.filter",
    )

    # Check the recorded comparison values
    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 1
    comp = session.query(db_orm.Comparison).one()
    query = "689d3fd6881db36b5e08329cf23cecdd"  # MGV-GENOME-0264574.fas
    subject = "78975d5144a1cd12e98898d573cf6536"  # MGV-GENOME-0266457.fna
    pytest.approx(
        comp.identity,
        get_matrix_entry(
            input_genomes_tiny / "matrices" / "ANIm_identity.tsv", query, subject
        ),
    )
    pytest.approx(
        comp.aln_length,
        get_matrix_entry(
            input_genomes_tiny / "matrices" / "ANIm_aln_lengths.tsv", query, subject
        ),
    )
    pytest.approx(
        comp.sim_errors,
        get_matrix_entry(
            input_genomes_tiny / "matrices" / "ANIm_sim_errors.tsv", query, subject
        ),
    )
    pytest.approx(
        comp.cov_query,
        get_matrix_entry(
            input_genomes_tiny / "matrices" / "ANIm_coverage.tsv", query, subject
        ),
    )
    session.close()
    tmp_db.unlink()


def test_logging_anim_bad_alignment(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_bad_alignments: Path,
) -> None:
    """Check can log a ANIm comparison to DB (bad alignments)."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_nucmer()

    private_cli.log_run(
        fasta=input_genomes_bad_alignments,
        database=tmp_db,
        cmdline="pyani-plus anim ...",
        status="Testing",
        name="Testing log_anim",
        method="ANIm",
        program=tool.exe_path.stem,
        version=tool.version,
        mode=method_anim.MODE,
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    private_cli.log_anim(
        database=tmp_db,
        run_id=1,
        query_fasta=input_genomes_bad_alignments / "MGV-GENOME-0264574.fas",
        subject_fasta=input_genomes_bad_alignments / "MGV-GENOME-0357962.fna",
        deltafilter=input_genomes_bad_alignments
        / "intermediates/ANIm/MGV-GENOME-0264574_vs_MGV-GENOME-0357962.filter",
    )

    # Check the recorded comparison values
    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 1
    comp = session.query(db_orm.Comparison).one()
    assert all(
        value is None
        for value in [
            comp.identity,
            comp.aln_length,
            comp.sim_errors,
            comp.cov_query,
            comp.cov_subject,
        ]
    )
    assert (
        comp.query_hash == "689d3fd6881db36b5e08329cf23cecdd"
    )  # MGV-GENOME-0264574.fas
    assert (
        comp.subject_hash == "a30481565b45f6bbc6ce5260503067e0"
    )  # MGV-GENOME-0357962.fna
    session.close()
    tmp_db.unlink()
