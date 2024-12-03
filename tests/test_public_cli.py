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

import filecmp
from pathlib import Path

import pandas as pd
import pytest

from pyani_plus import db_orm, public_cli, tools
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


def test_resume_empty(tmp_path: str) -> None:
    """Check list-runs with no data."""
    with pytest.raises(
        SystemExit, match="ERROR: Database /does/not/exist does not exist"
    ):
        public_cli.resume(database=Path("/does/not/exist"))

    tmp_db = Path(tmp_path) / "resume-empty.sqlite"
    session = db_orm.connect_to_db(tmp_db)
    session.close()

    with pytest.raises(
        SystemExit, match="ERROR: Database .*/resume-empty.sqlite contains no runs."
    ):
        public_cli.resume(database=tmp_db)

    with pytest.raises(
        SystemExit,
        match="ERROR: Database .*/resume-empty.sqlite has no run-id 1.",
    ):
        public_cli.resume(database=tmp_db, run_id=1)


def test_partial_run(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check list-runs and export-run with mock data including a partial run."""
    tmp_db = Path(tmp_path) / "list-runs.sqlite"
    session = db_orm.connect_to_db(tmp_db)
    config = db_orm.db_configuration(
        session, "fastANI", "fastani", "1.2.3", create=True
    )
    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(session, filename, md5, create=True)

    # Record 4 of the possible 9 comparisons:
    for query_hash in list(fasta_to_hash.values())[0:2]:
        for subject_hash in list(fasta_to_hash.values())[0:2]:
            db_orm.db_comparison(
                session,
                config.configuration_id,
                query_hash,
                subject_hash,
                1.0 if query_hash == subject_hash else 0.99,
                12345,
            )

    db_orm.add_run(
        session,
        config,
        cmdline="pyani fastani ...",
        fasta_directory=input_genomes_tiny,
        status="Empty",
        name="Trial A",
        fasta_to_hash={},
    )
    db_orm.add_run(
        session,
        config,
        cmdline="pyani fastani ...",
        fasta_directory=input_genomes_tiny,
        status="Partial",
        name="Trial B",
        fasta_to_hash=fasta_to_hash,  # 3/3 genomes, so only have 4/9 comparisons
    )
    db_orm.add_run(
        session,
        config,
        cmdline="pyani fastani ...",
        fasta_directory=input_genomes_tiny,
        status="Done",
        name="Trial C",
        # This run uses just 2/3 genomes, but has all 4/4 = 2*2 comparisons:
        fasta_to_hash=dict(list(fasta_to_hash.items())[0:2]),
    )

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 3 analysis runs in " in output, output
    assert " 0/0=0² │ Empty " in output, output
    assert " 4/9=3² │ Partial " in output, output
    assert " 4/4=2² │ Done " in output, output

    # Unlike a typical method calculation, we have not triggered
    # .cache_comparisons() yet, so that will happen in export_run.
    public_cli.export_run(database=tmp_db, run_id=3, outdir=tmp_path)
    output = capsys.readouterr().out
    assert f"Wrote matrices to {tmp_path}" in output, output

    # By construction run 2 is partial, only 4 of 9 matrix entries are
    # defined - this should fail
    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Database .*/list-runs.sqlite run-id 2 has only 4 of 3²=9 comparisons, 5 needed"
        ),
    ):
        public_cli.export_run(database=tmp_db, run_id=2, outdir=tmp_path)

    # Resuming the partial job should fail as the fastANI version won't match:
    with pytest.raises(
        SystemExit,
        match="ERROR: We have fastANI version .*, but run-id 2 used fastani version 1.2.3 instead.",
    ):
        public_cli.resume(database=tmp_db, run_id=2)

    # Resuming run 1 should fail as no genomes:
    with pytest.raises(
        SystemExit,
        match="ERROR: No genomes recorded for run-id 1, cannot resume.",
    ):
        public_cli.resume(database=tmp_db, run_id=1)

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
        session,
        config,
        cmdline="pyani fastani ...",
        fasta_directory=Path("/does/not/exist"),
        status="Empty",
        name="Trial A",
    )
    db_orm.add_run(
        session,
        config,
        cmdline="pyani fastani ...",
        fasta_directory=Path("/does/not/exist"),
        status="Empty",
        name="Trial B",
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


def test_anim(tmp_path: str, input_genomes_tiny: Path) -> None:
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
        input_genomes_tiny / "matrices" / "ANIm_identity.tsv",
        out / "ANIm_identity.tsv",
    )


def test_dnadiff(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check dnadiff run (default settings)."""
    out = Path(tmp_path)
    tmp_db = out / "example.sqlite"
    # Leaving out name, so can check the default worked
    public_cli.dnadiff(database=tmp_db, fasta=input_genomes_tiny, create_db=True)
    public_cli.export_run(database=tmp_db, outdir=out)
    # Fuzzy, 0.9963 from dnadiff tool != 0.9962661747 from our code
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "dnadiff_identity.tsv",
        out / "dnadiff_identity.tsv",
        atol=5e-5,
    )
    session = db_orm.connect_to_db(tmp_db)
    run = session.query(db_orm.Run).one()
    assert run.name == "3 genomes using dnadiff"
    session.close()


def test_anib(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check ANIb run (default settings)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "example.sqlite"
    public_cli.anib(
        database=tmp_db,
        fasta=input_genomes_tiny,
        name="Test Run",
        create_db=True,
        temp=tmp_dir,
    )

    for file in (input_genomes_tiny / "intermediates/ANIb").glob("*.f*"):
        assert filecmp.cmp(file, tmp_dir / file), f"Wrong fragmented FASTA {file.name}"

    # Could check the BLAST DB *.njs files here too...

    for file in (input_genomes_tiny / "intermediates/ANIb").glob("*_vs_*.tsv"):
        assert filecmp.cmp(file, tmp_dir / file), f"Wrong blastn output in {file.name}"

    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "ANIb_identity.tsv",
        tmp_dir / "ANIb_identity.tsv",
    )


def test_fastani(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check fastANI run (default settings)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "example.sqlite"
    public_cli.fastani(
        database=tmp_db,
        fasta=input_genomes_tiny,
        name="Test Run",
        create_db=True,
        temp=tmp_dir,
    )

    for file in (input_genomes_tiny / "intermediates/fastANI").glob("*_vs_*.fastani"):
        assert filecmp.cmp(file, tmp_dir / file), f"Wrong fastANI output in {file.name}"

    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "fastANI_identity.tsv",
        tmp_dir / "fastANI_identity.tsv",
    )


def test_sourmash(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check sourmash run (default settings except scaled=300)."""
    out = Path(tmp_path)
    tmp_db = out / "example.sqlite"
    public_cli.sourmash(
        database=tmp_db,
        fasta=input_genomes_tiny,
        name="Test Run",
        scaled=300,
        create_db=True,
    )
    public_cli.export_run(database=tmp_db, outdir=out)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "sourmash_identity.tsv",
        out / "sourmash_identity.tsv",
    )


def test_branchwater(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check sourmash run (default settings except scaled=300)."""
    out = Path(tmp_path)
    tmp_db = out / "example.sqlite"
    public_cli.branchwater(
        database=tmp_db,
        fasta=input_genomes_tiny,
        name="Test Run",
        scaled=300,
        create_db=True,
    )
    public_cli.export_run(database=tmp_db, outdir=out)
    # Should match the sourmash output (but computed quicker)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "sourmash_identity.tsv",
        out / "branchwater_identity.tsv",
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


def test_resume_partial_fastani(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check list-runs and export-run with mock data including a partial fastANI run."""
    tmp_db = Path(tmp_path) / "resume.sqlite"
    tool = tools.get_fastani()
    session = db_orm.connect_to_db(tmp_db)
    config = db_orm.db_configuration(
        session,
        "fastANI",
        tool.exe_path.stem,
        tool.version,
        kmersize=14,
        fragsize=2500,
        minmatch=0.3,
        create=True,
    )

    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(session, filename, md5, create=True)

    # Record 8 of the possible 9 comparisons:
    genomes = list(fasta_to_hash.values())
    for query_hash in genomes:
        for subject_hash in genomes:
            if query_hash == genomes[2] and subject_hash == genomes[0]:
                # Skip one comparison
                continue
            db_orm.db_comparison(
                session,
                config.configuration_id,
                query_hash,
                subject_hash,
                1.0 if query_hash is subject_hash else 0.99,
                12345,
            )

    db_orm.add_run(
        session,
        config,
        cmdline="pyani-plus fastani ...",
        fasta_directory=input_genomes_tiny,
        status="Partial",
        name="Test Resuming A Run",
        fasta_to_hash=fasta_to_hash,  # 3/3 genomes, so only have 8/9 comparisons
    )
    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " 8/9=3² │ Partial " in output, output

    public_cli.resume(database=tmp_db)
    output = capsys.readouterr().out
    assert "Resuming run-id 1, the only run" in output, output
    assert "Database already has 8 of 3²=9 comparisons, 1 needed" in output, output

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " 9/9=3² │ Done " in output, output


def test_resume_partial_anib(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check list-runs and export-run with mock data including a partial ANIb run."""
    tmp_db = Path(tmp_path) / "resume.sqlite"
    tool = tools.get_blastn()
    session = db_orm.connect_to_db(tmp_db)
    config = db_orm.db_configuration(
        session,
        "ANIb",
        tool.exe_path.stem,
        tool.version,
        fragsize=1234,
        create=True,
    )

    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(session, filename, md5, create=True)

    # Record 8 of the possible 9 comparisons:
    genomes = list(fasta_to_hash.values())
    for query_hash in genomes:
        for subject_hash in genomes:
            if query_hash == genomes[2] and subject_hash == genomes[0]:
                # Skip one comparison
                continue
            db_orm.db_comparison(
                session,
                config.configuration_id,
                query_hash,
                subject_hash,
                1.0 if query_hash is subject_hash else 0.99,
                12345,
            )

    db_orm.add_run(
        session,
        config,
        cmdline="pyani-plus anib ...",
        fasta_directory=input_genomes_tiny,
        status="Partial",
        name="Test Resuming A Run",
        fasta_to_hash=fasta_to_hash,  # 3/3 genomes, so only have 8/9 comparisons
    )
    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " 8/9=3² │ Partial " in output, output

    public_cli.resume(database=tmp_db)
    output = capsys.readouterr().out
    assert "Resuming run-id 1, the only run" in output, output
    assert "Database already has 8 of 3²=9 comparisons, 1 needed" in output, output

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " 9/9=3² │ Done " in output, output


def test_resume_partial_anim(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check list-runs and export-run with mock data including a partial ANIm run."""
    tmp_db = Path(tmp_path) / "resume.sqlite"
    tool = tools.get_nucmer()
    session = db_orm.connect_to_db(tmp_db)
    config = db_orm.db_configuration(
        session,
        "ANIm",
        tool.exe_path.stem,
        tool.version,
        mode="maxmatch",
        create=True,
    )

    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(session, filename, md5, create=True)

    # Record 8 of the possible 9 comparisons:
    genomes = list(fasta_to_hash.values())
    for query_hash in genomes:
        for subject_hash in genomes:
            if query_hash == genomes[2] and subject_hash == genomes[0]:
                # Skip one comparison
                continue
            db_orm.db_comparison(
                session,
                config.configuration_id,
                query_hash,
                subject_hash,
                1.0 if query_hash is subject_hash else 0.99,
                12345,
            )

    db_orm.add_run(
        session,
        config,
        cmdline="pyani-plus ANIm --mode maxmatch ...",
        fasta_directory=input_genomes_tiny,
        status="Partial",
        name="Test Resuming A Run",
        fasta_to_hash=fasta_to_hash,  # 3/3 genomes, so only have 8/9 comparisons
    )
    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " 8/9=3² │ Partial " in output, output

    public_cli.resume(database=tmp_db)
    output = capsys.readouterr().out
    assert "Resuming run-id 1, the only run" in output, output
    assert "Database already has 8 of 3²=9 comparisons, 1 needed" in output, output

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " 9/9=3² │ Done " in output, output


def test_resume_partial_sourmash(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check list-runs and export-run with mock data including a partial sourmash run."""
    tmp_db = Path(tmp_path) / "resume.sqlite"
    tool = tools.get_sourmash()
    session = db_orm.connect_to_db(tmp_db)
    config = db_orm.db_configuration(
        session,
        "sourmash",
        tool.exe_path.stem,
        tool.version,
        kmersize=31,  # must be 31 to match the sig files in fixtures
        extra="scaled=300",
        mode="containment",
        create=True,
    )

    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(session, filename, md5, create=True)

    # Record 4 of the possible 9 comparisons,
    # mimicking what might happen when a 2x2 run is expanded to 3x3
    genomes = list(fasta_to_hash.values())
    for query_hash in genomes[:-1]:
        for subject_hash in genomes[:-1]:
            db_orm.db_comparison(
                session,
                config.configuration_id,
                query_hash,
                subject_hash,
                1.0 if query_hash is subject_hash else 0.99,
            )

    db_orm.add_run(
        session,
        config,
        cmdline="pyani-plus sourmash ...",
        fasta_directory=input_genomes_tiny,
        status="Partial",
        name="Test Resuming A Run",
        fasta_to_hash=fasta_to_hash,  # all 3/3 genomes, but only have 4/9 comparisons
    )
    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " 4/9=3² │ Partial " in output, output

    public_cli.resume(database=tmp_db)
    output = capsys.readouterr().out
    assert "Resuming run-id 1, the only run" in output, output
    assert "Database already has 4 of 3²=9 comparisons, 5 needed" in output, output

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " 9/9=3² │ Done " in output, output


def test_resume_partial_branchwater(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check list-runs and export-run with mock data including a partial sourmash run."""
    tmp_db = Path(tmp_path) / "resume.sqlite"
    tool = tools.get_sourmash()
    session = db_orm.connect_to_db(tmp_db)
    config = db_orm.db_configuration(
        session,
        "branchwater",
        tool.exe_path.stem,
        tool.version,
        kmersize=31,  # must be 31 to match the sig files in fixtures
        extra="scaled=300",
        mode="containment",
        create=True,
    )

    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(session, filename, md5, create=True)

    # Record 4 of the possible 9 comparisons,
    # mimicking what might happen when a 2x2 run is expanded to 3x3
    genomes = list(fasta_to_hash.values())
    for query_hash in genomes[:-1]:
        for subject_hash in genomes[:-1]:
            db_orm.db_comparison(
                session,
                config.configuration_id,
                query_hash,
                subject_hash,
                1.0 if query_hash is subject_hash else 0.99,
            )

    db_orm.add_run(
        session,
        config,
        cmdline="pyani-plus branchwater ...",
        fasta_directory=input_genomes_tiny,
        status="Partial",
        name="Test Resuming A Run",
        fasta_to_hash=fasta_to_hash,  # all 3/3 genomes, but only have 4/9 comparisons
    )
    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " 4/9=3² │ Partial " in output, output

    public_cli.resume(database=tmp_db)
    output = capsys.readouterr().out
    assert "Resuming run-id 1, the only run" in output, output
    assert "Database already has 4 of 3²=9 comparisons, 5 needed" in output, output

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " 9/9=3² │ Done " in output, output


def test_resume_dir_gone(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check expected failure trying to resume without the input directory."""
    tmp_db = Path(tmp_path) / "resume.sqlite"
    tool = tools.get_fastani()
    session = db_orm.connect_to_db(tmp_db)
    config = db_orm.db_configuration(
        session,
        "fastANI",
        tool.exe_path.stem,
        tool.version,
        kmersize=14,
        fragsize=2500,
        minmatch=0.3,
        create=True,
    )

    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(session, filename, md5, create=True)

    db_orm.add_run(
        session,
        config,
        cmdline="pyani ...",
        fasta_directory=Path("/mnt/shared/old"),
        status="Partial",
        name="Test how resuming without input directory present fails",
        fasta_to_hash=fasta_to_hash,
    )
    session.close()

    with pytest.raises(
        SystemExit,
        match=r"ERROR: run-id 1 used input folder /mnt/shared/old, but that is not a directory \(now\).",
    ):
        public_cli.resume(database=tmp_db)


def test_resume_unknown(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check expected failure trying to resume an unknown method."""
    tmp_db = Path(tmp_path) / "resume.sqlite"
    session = db_orm.connect_to_db(tmp_db)
    config = db_orm.db_configuration(
        session,
        "guessing",
        "gestimator",
        "1.2.3b4",
        kmersize=51,
        fragsize=999,
        minmatch=0.25,
        create=True,
    )

    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(session, filename, md5, create=True)

    db_orm.add_run(
        session,
        config,
        cmdline="pyani ...",
        fasta_directory=input_genomes_tiny,
        status="Partial",
        name="Test how resuming from an unknown method fails",
        fasta_to_hash=fasta_to_hash,
    )
    session.close()

    with pytest.raises(
        SystemExit,
        match="ERROR: Unknown method guessing for run-id 1 in .*/resume.sqlite",
    ):
        public_cli.resume(database=tmp_db)


def test_resume_complete(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check resume works for all the methods (using completed runs for speed)."""
    tmp_db = Path(tmp_path) / "resume.sqlite"
    session = db_orm.connect_to_db(tmp_db)

    for index, (method, tool) in enumerate(
        [
            ("ANIm", tools.get_nucmer()),
            ("dnadiff", tools.get_nucmer()),
            ("ANIb", tools.get_blastn()),
            ("fastANI", tools.get_fastani()),
            ("sourmash", tools.get_sourmash()),
        ]
    ):
        config = db_orm.db_configuration(
            session,
            method,
            tool.exe_path.stem,
            tool.version,
            create=True,
        )

        fasta_to_hash = {
            filename: file_md5sum(filename)
            for filename in sorted(input_genomes_tiny.glob("*.f*"))
        }
        for filename, md5 in fasta_to_hash.items():
            db_orm.db_genome(session, filename, md5, create=True)

        # Record dummy values for all of the possible 9 comparisons:
        for query_hash in fasta_to_hash.values():
            for subject_hash in fasta_to_hash.values():
                db_orm.db_comparison(
                    session,
                    config.configuration_id,
                    query_hash,
                    subject_hash,
                    1.0 if query_hash == subject_hash else 0.99,
                    12345,
                )

        db_orm.add_run(
            session,
            config,
            cmdline=f"pyani {method} --database ...",
            fasta_directory=input_genomes_tiny,
            status="Complete",
            name=f"Test resuming a complete {method} run",
            fasta_to_hash=fasta_to_hash,
        )
        public_cli.resume(database=tmp_db)
        output = capsys.readouterr().out
        if index == 0:
            assert "Resuming run-id 1, the only run" in output, output
        else:
            assert f"Resuming run-id {index+1}, the latest run" in output, output
        assert "Database already has all 3²=9 comparisons" in output, output


def test_resume_fasta_gone(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check resume error handling when a FASTA file is missing."""
    tmp_indir = Path(tmp_path) / "input"
    tmp_indir.mkdir()
    tmp_db = Path(tmp_path) / "resume.sqlite"
    session = db_orm.connect_to_db(tmp_db)
    tool = tools.get_blastn()
    config = db_orm.db_configuration(
        session,
        "ANIb",
        tool.exe_path.stem,
        tool.version,
        fragsize=1500,
        create=True,
    )

    # Record the genome entries under input_genomes_tiny
    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(session, filename, md5, create=True)

    # Setup a subset under tmp_dir
    for filename in list(fasta_to_hash)[:-1]:
        fasta = Path(filename)
        (tmp_indir / fasta.name).symlink_to(fasta)
    assert len(fasta_to_hash) - 1 == len(list(tmp_indir.glob("*.f*"))), list(
        tmp_indir.glob("*.f*")
    )
    missing = list(fasta_to_hash)[-1]

    # Record dummy values for all of the possible 9 comparisons:
    for query_hash in fasta_to_hash.values():
        for subject_hash in fasta_to_hash.values():
            db_orm.db_comparison(
                session,
                config.configuration_id,
                query_hash,
                subject_hash,
                1.0 if query_hash == subject_hash else 0.99,
                12345,
            )

    db_orm.add_run(
        session,
        config,
        cmdline="pyani ANIb --database ...",
        fasta_directory=tmp_indir,  # NOT using input_genomes_tiny
        status="Complete",
        name="Test resuming when the FASTA files have gone",
        fasta_to_hash=fasta_to_hash,
    )

    # Won't work with a FASTA missing:
    with pytest.raises(
        SystemExit,
        match=(
            f"ERROR: run-id 1 used .*/{Path(missing).name} with MD5 {fasta_to_hash[missing]}"
            " but this FASTA file no longer exists"
        ),
    ):
        public_cli.resume(database=tmp_db)

    # Should work with all the FASTA files, even though now in a different directory
    # to that logged in the genome table (as could happen via an older run):
    (tmp_indir / Path(missing).name).symlink_to(Path(missing))
    public_cli.resume(database=tmp_db)
    output = capsys.readouterr().out
    assert "Database already has all 3²=9 comparisons" in output, output

    # Should work even with extra FASTA files, real world use case:
    # Using input directory like /mnt/shared/genomes
    # Start a run when there are 100 genomes (but it was not completed)
    # Add some more genomes, say now 150 genomes
    # Resume the run - it should only operate on the first 100 genomes!
    with (tmp_indir / "extra.fasta").open("w") as handle:
        handle.write(">recently-added-genome\nACGTACGTAGT\n")
    public_cli.resume(database=tmp_db)
    output = capsys.readouterr().out
    assert "Database already has all 3²=9 comparisons" in output, output
