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
"""Tests for the pyani_plus/db_orm.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pathlib import Path

import pandas as pd
import pytest

from pyani_plus import db_orm, public_cli
from pyani_plus.utils import file_md5sum


# This is very similar to the functions under tests/snakemake/__init__.py
def compare_matrix_files(
    expected_file: Path, new_file: Path, atol: float | None = None
) -> None:
    """Compare two matrix files (after sorting)."""
    assert expected_file.is_file(), f"Missing expected {expected_file}"
    assert new_file.is_file(), f"Missing output {new_file}"
    expected_df = (
        pd.read_csv(expected_file, sep="\t", header=0, index_col=0)
        .sort_index(axis=0)
        .sort_index(axis=1)
    )
    new_df = (
        pd.read_csv(new_file, sep="\t", header=0, index_col=0)
        .sort_index(axis=0)
        .sort_index(axis=1)
    )
    if atol is None:
        pd.testing.assert_frame_equal(expected_df, new_df, obj=new_file)
    else:
        pd.testing.assert_frame_equal(expected_df, new_df, obj=new_file, atol=atol)


def test_check_db() -> None:
    """Check check_db error conditions."""
    with pytest.raises(
        SystemExit,
        match="ERROR: Database /does/not/exist does not exist, but not using --create-db",
    ):
        public_cli.check_db(Path("/does/not/exist"), create_db=False)

    # This is fine:
    public_cli.check_db(":memory:", create_db=True)


def test_check_fasta(tmp_path: str) -> None:
    """Check error conditions."""
    with pytest.raises(
        SystemExit, match="ERROR: FASTA input /does/not/exist is not a directory"
    ):
        public_cli.check_fasta(Path("/does/not/exist"))

    with pytest.raises(SystemExit, match="ERROR: No FASTA input genomes under "):
        public_cli.check_fasta(Path(tmp_path))


def test_list_runs_empty(capsys: pytest.CaptureFixture[str], tmp_path: str) -> None:
    """Check list-runs with no data."""
    with pytest.raises(
        SystemExit, match="ERROR: Database /does/not/exist does not exist"
    ):
        public_cli.list_runs(database=Path("/does/not/exist"))

    tmp_db = Path(tmp_path) / "list-runs-empty.sqlite"
    session = db_orm.connect_to_db(tmp_db)
    session.close()

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 0 analysis runs in " in output, output


def test_list_runs(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check list-runs with mock data."""
    tmp_db = Path(tmp_path) / "list-runs.sqlite"
    session = db_orm.connect_to_db(tmp_db)
    config = db_orm.db_configuration(
        session, "fastANI", "fastani", "1.2.3", create=True
    )
    genomes = [
        db_orm.db_genome(session, filename, file_md5sum(filename), create=True)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    ]
    # Record 4 of the possible 9 comparisons:
    for query in genomes[0:2]:
        for subject in genomes[0:2]:
            db_orm.db_comparison(
                session,
                config.configuration_id,
                query.genome_hash,
                subject.genome_hash,
                1.0 if query is subject else 0.99,
                12345,
            )
    db_orm.add_run(
        session,
        config,
        cmdline="pyani fastani ...",
        status="Empty",
        name="Trial A",
        genomes=[],
    )
    db_orm.add_run(
        session,
        config,
        cmdline="pyani fastani ...",
        status="Partial",
        name="Trial B",
        genomes=genomes,  # 3/3 genomes, so only have 4/9 comparisons
    )
    db_orm.add_run(
        session,
        config,
        cmdline="pyani fastani ...",
        status="Done",
        name="Trial C",
        genomes=genomes[0:2],  # 2/3 genomes, but have all 4/4 comparisons
    )
    # Can we test stdout, should say 3 runs:
    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 3 analysis runs in " in output, output
    assert " 0/0=0² │ Empty " in output, output
    assert " 4/9=3² │ Partial " in output, output
    assert " 4/4=2² │ Done " in output, output
    tmp_db.unlink()


def test_export_run(tmp_path: str) -> None:
    """Check export run failures."""
    tmp = Path(tmp_path)

    with pytest.raises(
        SystemExit, match="ERROR: Database /does/not/exist does not exist"
    ):
        public_cli.export_run(database=Path("/does/not/exist"), outdir=tmp)

    with pytest.raises(
        SystemExit, match="ERROR: Output directory /does/not/exist does not exist"
    ):
        public_cli.export_run(database=":memory:", outdir=Path("/does/not/exist/"))

    tmp_db = Path(tmp_path) / "export.sqlite"
    session = db_orm.connect_to_db(tmp_db)
    with pytest.raises(
        SystemExit, match="ERROR: Database .*/export.sqlite contains no runs."
    ):
        public_cli.export_run(database=tmp_db, outdir=tmp)

    config = db_orm.db_configuration(
        session, "fastANI", "fastani", "1.2.3", create=True
    )
    db_orm.add_run(
        session, config, cmdline="pyani fastani ...", status="Empty", name="Trial A"
    )
    db_orm.add_run(
        session, config, cmdline="pyani fastani ...", status="Empty", name="Trial B"
    )
    with pytest.raises(
        SystemExit,
        match="ERROR: Database .*/export.sqlite contains 2 runs, use --run-id to specify which.",
    ):
        public_cli.export_run(database=tmp_db, outdir=tmp)
    with pytest.raises(
        SystemExit, match="ERROR: Database .*/export.sqlite has no run-id 3."
    ):
        public_cli.export_run(database=tmp_db, outdir=tmp, run_id=3)
    with pytest.raises(
        SystemExit,
        match="ERROR: Database .*/export.sqlite run-id 1 has no comparisons",
    ):
        public_cli.export_run(database=tmp_db, outdir=tmp, run_id=1)
    tmp_db.unlink()


def test_anim(tmp_path: str, input_genomes_tiny: Path, dir_anim_results: Path) -> None:
    """Check ANIm run."""
    out = Path(tmp_path)
    tmp_db = out / "example.sqlite"
    public_cli.anim(
        database=tmp_db, fasta=input_genomes_tiny, name="Test Run", create_db=True
    )
    # Now do it again - it should reuse the calculations:
    public_cli.anim(
        database=tmp_db, fasta=input_genomes_tiny, name="Test Run", create_db=False
    )
    public_cli.export_run(database=tmp_db, outdir=out, run_id=1)  # have two runs
    compare_matrix_files(
        dir_anim_results / "matrix_identity.tsv", out / "ANIm_identity.tsv"
    )


def test_dnadiff(
    tmp_path: str, input_genomes_tiny: Path, dir_dnadiff_matrices: Path
) -> None:
    """Check dnadiff run (default settings)."""
    out = Path(tmp_path)
    tmp_db = out / "example.sqlite"
    public_cli.dnadiff(
        database=tmp_db, fasta=input_genomes_tiny, name="Test Run", create_db=True
    )
    public_cli.export_run(database=tmp_db, outdir=out)
    # Fuzzy, 0.9963 from dnadiff tool != 0.9962661747 from our code
    compare_matrix_files(
        dir_dnadiff_matrices / "matrix_identity.tsv",
        out / "dnadiff_identity.tsv",
        atol=5e-5,
    )


def test_anib(tmp_path: str, input_genomes_tiny: Path, dir_anib_results: Path) -> None:
    """Check ANIb run (default settings)."""
    out = Path(tmp_path)
    tmp_db = out / "example.sqlite"
    public_cli.anib(
        database=tmp_db, fasta=input_genomes_tiny, name="Test Run", create_db=True
    )
    public_cli.export_run(database=tmp_db, outdir=out)
    compare_matrix_files(
        dir_anib_results / "matrix_identity.tsv", out / "ANIb_identity.tsv"
    )


def test_fastani(
    tmp_path: str, input_genomes_tiny: Path, fastani_matrices: Path
) -> None:
    """Check fastANI run (default settings)."""
    out = Path(tmp_path)
    tmp_db = out / "example.sqlite"
    public_cli.fastani(
        database=tmp_db, fasta=input_genomes_tiny, name="Test Run", create_db=True
    )
    public_cli.export_run(database=tmp_db, outdir=out)
    compare_matrix_files(
        fastani_matrices / "matrix_identity.tsv", out / "fastANI_identity.tsv"
    )


def test_fastani_dups(tmp_path: str) -> None:
    """Check fastANI run (duplicate FASTA inputs)."""
    tmp = Path(tmp_path)
    tmp_db = tmp / "example.sqlite"
    for name in ("alpha", "beta", "gamma"):
        with (tmp / (name + ".fasta")).open("w") as handle:
            handle.write(">genome\nACGTACGT\n")
    with pytest.raises(
        SystemExit, match="ERROR - Multiple genomes with same MD5 checksum"
    ):
        public_cli.fastani(
            database=tmp_db, fasta=tmp, name="Test duplicates fail", create_db=True
        )
