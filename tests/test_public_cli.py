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
"""Tests for the pyani_plus/public_cli.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

import filecmp
import gzip
import logging
import shutil
from pathlib import Path

import pandas as pd
import pytest

from pyani_plus import GRAPHICS_FORMATS, db_orm, public_cli, setup_logger, tools
from pyani_plus.public_cli_args import ToolExecutor
from pyani_plus.utils import file_md5sum


@pytest.fixture(scope="session")
def gzipped_tiny_example(
    tmp_path_factory: pytest.TempPathFactory, input_genomes_tiny: Path
) -> Path:
    """Make a gzipped version of the viral directory of FASTA files."""
    gzip_dir = tmp_path_factory.mktemp("gzipped_" + input_genomes_tiny.stem)
    for fasta in input_genomes_tiny.glob("*.f*"):
        with (
            (input_genomes_tiny / fasta).open("rb") as f_in,
            gzip.open(str(gzip_dir / (fasta.name + ".gz")), "wb") as f_out,
        ):
            shutil.copyfileobj(f_in, f_out)
    return gzip_dir


@pytest.fixture(scope="session")
def evil_example(
    tmp_path_factory: pytest.TempPathFactory, input_genomes_tiny: Path
) -> Path:
    """Make a version of the viral directory of FASTA files using spaces, emoji, etc."""
    space_dir = tmp_path_factory.mktemp("with spaces " + input_genomes_tiny.stem)
    for fasta in input_genomes_tiny.glob("*.f*"):
        space_fasta = space_dir / ("🦠 : " + fasta.name.replace("-", " "))
        space_fasta.symlink_to(fasta)
    return space_dir


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
    logger = setup_logger(None)
    with pytest.raises(
        SystemExit,
        match="Database /does/not/exist does not exist, but not using --create-db",
    ):
        public_cli.check_db(logger, Path("/does/not/exist"), create_db=False)

    # This is fine:
    public_cli.check_db(logger, ":memory:", create_db=True)


def test_check_fasta(tmp_path: str) -> None:
    """Check error conditions."""
    tmp_dir = Path(tmp_path)
    logger = setup_logger(Path("-"))
    with pytest.raises(
        SystemExit, match="FASTA input /does/not/exist is not a directory"
    ):
        public_cli.check_fasta(logger, Path("/does/not/exist"))

    with pytest.raises(SystemExit, match="No FASTA input genomes under "):
        public_cli.check_fasta(logger, tmp_dir)


def test_check_start_and_run(tmp_path: str) -> None:
    """Check error conditions."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "x.db"
    (tmp_dir / "broken.fasta").symlink_to("/does/not/exist/example.fna")
    logger = setup_logger(Path("-"))
    with pytest.raises(SystemExit, match="Input /.*/broken.fasta is a broken symlink"):
        public_cli.start_and_run_method(
            logger,
            ToolExecutor.local,
            Path(),
            tmp_dir,
            None,
            tmp_db,
            tmp_dir / "x.log",
            "test",
            "guessing",
            tmp_dir,
            None,
        )


def test_delete_empty(tmp_path: str) -> None:
    """Check delete-run with no data."""
    with pytest.raises(SystemExit, match="Database /does/not/exist does not exist"):
        public_cli.delete_run(database=Path("/does/not/exist"), log=Path("-"))

    tmp_db = Path(tmp_path) / "list-runs-empty.sqlite"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    session.close()

    with pytest.raises(SystemExit, match="Database contains no runs."):
        public_cli.delete_run(database=tmp_db)

    with pytest.raises(SystemExit, match="Database has no run-id 1."):
        public_cli.resume(database=tmp_db, run_id=1)


def test_list_runs_empty(capsys: pytest.CaptureFixture[str], tmp_path: str) -> None:
    """Check list-runs with no data."""
    with pytest.raises(SystemExit, match="Database /does/not/exist does not exist"):
        public_cli.list_runs(database=Path("/does/not/exist"))

    tmp_db = Path(tmp_path) / "list runs empty.sqlite"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    session.close()

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 0 analysis runs in " in output, output


def test_resume_empty(tmp_path: str) -> None:
    """Check list-runs with no data."""
    with pytest.raises(SystemExit, match="Database /does/not/exist does not exist"):
        public_cli.resume(database=Path("/does/not/exist"))

    tmp_db = Path(tmp_path) / "resume-empty.sqlite"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    session.close()

    with pytest.raises(SystemExit, match="Database contains no runs."):
        public_cli.resume(database=tmp_db)

    with pytest.raises(SystemExit, match="Database has no run-id 1."):
        public_cli.resume(database=tmp_db, run_id=1)


def test_partial_run(  # noqa: PLR0915
    caplog: pytest.LogCaptureFixture,
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check list-runs and export-run with mock data including a partial run."""
    caplog.set_level(logging.INFO)
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "list runs.sqlite"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    config = db_orm.db_configuration(
        session, "fastANI", "fastani", "1.2.3", create=True
    )
    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(logger, session, filename, md5, create=True)

    # Record 4 of the possible 9 comparisons:
    for query_hash in list(fasta_to_hash.values())[1:]:
        for subject_hash in list(fasta_to_hash.values())[1:]:
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
        status="Running",  # simulate a partial run not ended cleanly
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
        fasta_to_hash=dict(list(fasta_to_hash.items())[1:]),
    )

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 3 analysis runs in " in output, output
    assert " Method  ┃ Done ┃ Null ┃ Miss ┃ Total ┃ Status " in output, output
    assert " fastANI │    0 │    0 │    0 │  0=0² │ Empty " in output, output
    assert " fastANI │    4 │    0 │    5 │  9=3² │ Running " in output, output
    assert " fastANI │    4 │    0 │    0 │  4=2² │ Done " in output, output

    # Unlike a typical method calculation, we have not triggered
    # .cache_comparisons() yet, so that will happen in export_run.
    caplog.clear()
    public_cli.export_run(database=tmp_db, run_id=3, outdir=tmp_path, label="md5")
    output = caplog.text
    assert f"Wrote matrices to {tmp_path}" in output, output
    with (tmp_dir / "fastANI_identity.tsv").open() as handle:
        assert (
            handle.readline()
            == "\t5584c7029328dc48d33f95f0a78f7e57\t78975d5144a1cd12e98898d573cf6536\n"
        )

    caplog.clear()
    public_cli.export_run(database=tmp_db, run_id=3, outdir=tmp_path, label="stem")
    output = caplog.text
    assert f"Wrote matrices to {tmp_path}" in output, output

    with (tmp_dir / "fastANI_identity.tsv").open() as handle:
        assert handle.readline() == "\tMGV-GENOME-0266457\tOP073605\n"

    caplog.clear()
    public_cli.export_run(database=tmp_db, run_id=3, outdir=tmp_path, label="filename")
    output = caplog.text
    assert f"Wrote matrices to {tmp_path}" in output, output
    with (tmp_dir / "fastANI_identity.tsv").open() as handle:
        assert handle.readline() == "\tMGV-GENOME-0266457.fna\tOP073605.fasta\n"
    caplog.clear()  # clear the above commands

    # By construction run 2 is partial, only 4 of 9 matrix entries are
    # defined - this should fail
    with pytest.raises(
        SystemExit,
        match=("run-id 2 has only 4 of 3²=9 comparisons, 5 needed"),
    ):
        public_cli.export_run(database=tmp_db, run_id=2, outdir=tmp_path)
    long_file = tmp_dir / "fastANI_run_2.tsv"
    assert long_file.is_file()
    with long_file.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        assert header == [
            "#Query",
            "Subject",
            "Identity",
            "Query-Cov",
            "Subject-Cov",
            "Hadamard",
            "tANI",
            "Align-Len",
            "Sim-Errors",
        ]
        assert sum(1 for _ in handle) == 4  # noqa: PLR2004

    # Resuming the partial job should fail as the fastANI version won't match:
    with pytest.raises(
        SystemExit,
        match="We have fastANI version .*, but run-id 2 used fastani version 1.2.3 instead.",
    ):
        public_cli.resume(database=tmp_db, run_id=2)

    # Resuming run 1 should fail as no genomes:
    with pytest.raises(
        SystemExit,
        match="No genomes recorded for run-id 1, cannot resume.",
    ):
        public_cli.resume(database=tmp_db, run_id=1)

    output = capsys.readouterr().out

    # Now delete the runs, and confirm what is left behind...
    caplog.clear()
    public_cli.delete_run(database=tmp_db, run_id=1)
    output = caplog.text
    assert "Run 1 contains 0/0=0² fastANI comparisons, status: Empty\n" in output

    # Forcing as this has data
    caplog.clear()
    public_cli.delete_run(database=tmp_db, run_id=3, force=True)
    output = caplog.text
    assert "Run 3 contains all 4=2² fastANI comparisons, status: Done\n" in output

    # Finally delete run 2, this will assume last one remaining
    caplog.clear()
    public_cli.delete_run(database=tmp_db, force=True)
    output = caplog.text
    assert "Deleting most recent run" in output, output
    assert "Run 2 contains 4/9=3² fastANI comparisons, status: Running" in output, (
        output
    )
    assert "Run name: Trial B" in output, output
    assert "Deleting a run still being computed will cause it to fail!" in output, (
        output
    )

    assert session.query(db_orm.Run).count() == 0
    assert session.query(db_orm.RunGenomeAssociation).count() == 0
    assert session.query(db_orm.Configuration).count() == 1
    assert session.query(db_orm.Genome).count() == 3  # noqa: PLR2004
    assert session.query(db_orm.Comparison).count() == 4  # noqa: PLR2004
    session.close()
    tmp_db.unlink()


def test_export_run_failures(tmp_path: str) -> None:
    """Check export run failures."""
    tmp_dir = Path(tmp_path)

    with pytest.raises(SystemExit, match="Database /does/not/exist does not exist"):
        public_cli.export_run(database=Path("/does/not/exist"), outdir=tmp_dir)

    tmp_db = tmp_dir / "empty.sqlite"
    tmp_db.touch()
    with pytest.raises(SystemExit, match="Database contains no runs."):
        public_cli.export_run(database=tmp_db, outdir=tmp_dir)

    tmp_db = tmp_dir / "export.sqlite"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
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
    with pytest.raises(SystemExit, match="Database has no run-id 3."):
        public_cli.export_run(database=tmp_db, outdir=tmp_dir, run_id=3)
    with pytest.raises(
        SystemExit,
        match="run-id 1 has no comparisons",
    ):
        public_cli.export_run(database=tmp_db, outdir=tmp_dir, run_id=1)
    # Should default to latest run, run-id 2
    with pytest.raises(
        SystemExit,
        match="run-id 2 has no comparisons",
    ):
        public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    tmp_db.unlink()


def test_export_duplicate_stem(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check export and plot run with duplicated stems.

    This should not happen naturally, it will fail via public CLI.
    """
    tmp_dir = Path(tmp_path)
    tmp_fasta = tmp_dir / "genomes"
    tmp_fasta.mkdir()
    (tmp_fasta / "example.fasta").symlink_to(
        input_genomes_tiny / "OP073605.fasta",
    )
    (tmp_fasta / "example.fna").symlink_to(
        input_genomes_tiny / "MGV-GENOME-0266457.fna",
    )
    (tmp_fasta / "example.fas").symlink_to(
        input_genomes_tiny / "MGV-GENOME-0264574.fas",
    )

    tmp_db = tmp_dir / "dup-stems.db"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    config = db_orm.db_configuration(
        session, "fastANI", "fastani", "1.2.3", create=True
    )
    fasta_to_hash = {fasta: file_md5sum(fasta) for fasta in tmp_fasta.glob("*.fa*")}
    for fasta, md5 in fasta_to_hash.items():
        db_orm.db_genome(logger, session, fasta, md5, create=True)
    db_orm.add_run(
        session,
        config,
        cmdline="pyani fastani ...",
        fasta_directory=input_genomes_tiny,
        status="Done",
        name="Trial B",
        fasta_to_hash=fasta_to_hash,
    )
    for query_hash in list(fasta_to_hash.values()):
        for subject_hash in list(fasta_to_hash.values()):
            db_orm.db_comparison(
                session,
                config.configuration_id,
                query_hash,
                subject_hash,
                1.0 if query_hash == subject_hash else 0.99,
                12345,
                cov_subject=1.0 if query_hash == subject_hash else 0.95,
                cov_query=1.0 if query_hash == subject_hash else 0.95,
            )
    session.commit()
    session.close()

    with pytest.raises(
        SystemExit,
        match="Duplicate filename stems, consider using MD5 labelling.",
    ):
        public_cli.export_run(database=tmp_db, outdir=tmp_dir / "out1")

    with pytest.raises(
        SystemExit,
        match="Duplicate filename stems, consider using MD5 labelling.",
    ):
        public_cli.plot_run(database=tmp_db, outdir=tmp_dir / "out2")
    with pytest.raises(
        SystemExit,
        match="Duplicate filename stems, consider using MD5 labelling.",
    ):
        public_cli.cli_classify(database=tmp_db, outdir=tmp_dir / "out3")


def test_plot_run_failures(tmp_path: str) -> None:
    """Check plot run failures."""
    tmp_dir = Path(tmp_path)

    with pytest.raises(SystemExit, match="Database /does/not/exist does not exist"):
        public_cli.plot_run(database=Path("/does/not/exist"), outdir=tmp_dir)

    logger = setup_logger(None)
    tmp_db = tmp_dir / "export.sqlite"
    session = db_orm.connect_to_db(logger, tmp_db)
    with pytest.raises(SystemExit, match="Database contains no runs."):
        public_cli.plot_run(database=tmp_db, outdir=tmp_dir)

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
    with pytest.raises(SystemExit, match="Database has no run-id 3."):
        public_cli.plot_run(database=tmp_db, outdir=tmp_dir, run_id=3)
    with pytest.raises(
        SystemExit,
        match="run-id 1 has no comparisons",
    ):
        public_cli.plot_run(database=tmp_db, outdir=tmp_dir, run_id=1)
    # Should default to latest run, run-id 2
    with pytest.raises(
        SystemExit,
        match="run-id 2 has no comparisons",
    ):
        public_cli.plot_run(database=tmp_db, outdir=tmp_dir)


def test_plot_run_comp_failures(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check plot-run-comp failures."""
    tmp_dir = Path(tmp_path)

    with pytest.raises(SystemExit, match="Database /does/not/exist does not exist"):
        public_cli.plot_run_comp(
            database=Path("/does/not/exist"), outdir=tmp_dir, run_ids="1,2"
        )

    tmp_db = tmp_dir / "export.sqlite"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    with pytest.raises(
        SystemExit,
        match="Database has no run-id 1. Use the list-runs command for more information.",
    ):
        public_cli.plot_run_comp(database=tmp_db, outdir=tmp_dir, run_ids="1,2")

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
    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(logger, session, filename, md5, create=True)
    genomes = list(fasta_to_hash.values())
    for query_hash in genomes:
        for subject_hash in genomes:
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
        cmdline="pyani fastani ...",
        fasta_directory=Path("/does/not/exist"),
        status="Empty",
        name="Trial B",
        fasta_to_hash=fasta_to_hash,
    )

    with pytest.raises(SystemExit, match="Run 1 has no comparisons"):
        public_cli.plot_run_comp(database=tmp_db, outdir=tmp_dir, run_ids="1,2")

    with pytest.raises(SystemExit, match="Need at least two runs for a comparison"):
        public_cli.plot_run_comp(database=tmp_db, outdir=tmp_dir, run_ids="2")

    with pytest.raises(
        SystemExit, match="Expected comma separated list of runs, not: 2-1"
    ):
        public_cli.plot_run_comp(database=tmp_db, outdir=tmp_dir, run_ids="2-1")

    with pytest.raises(SystemExit, match="Runs 2 and 1 have no comparisons in common"):
        public_cli.plot_run_comp(database=tmp_db, outdir=tmp_dir, run_ids="2,1")


def test_anim(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
    evil_example: Path,
) -> None:
    """Check ANIm run."""
    caplog.set_level(logging.INFO)
    tmp_dir = Path(tmp_path)
    # DB name with spaces, single quotes, emoji, in path & filename
    (tmp_dir / "user's 🔎 output").mkdir()
    tmp_db = tmp_dir / "user's 🔎 output" / "anim's 📦.sqlite"
    public_cli.cli_anim(
        database=tmp_db,
        fasta=evil_example,
        name="Spaces etc",
        create_db=True,
    )
    output = caplog.text
    assert "Database already has 0 of 3²=9 ANIm comparisons, 9 needed\n" in output
    logger = setup_logger(None)

    session = db_orm.connect_to_db(logger, tmp_db)
    count = session.query(db_orm.Comparison).count()
    session.close()

    # Now do it again - it should reuse the calculations:
    caplog.clear()
    public_cli.cli_anim(
        database=tmp_db,
        fasta=input_genomes_tiny,
        name="Simple names",
        create_db=False,
    )
    output = caplog.text
    assert "Database already has all 3²=9 ANIm comparisons\n" in output

    session = db_orm.connect_to_db(logger, tmp_db)
    assert count == session.query(db_orm.Comparison).count()
    session.close()

    public_cli.export_run(database=tmp_db, outdir=tmp_dir, run_id=2)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "ANIm_identity.tsv",
        tmp_dir / "ANIm_identity.tsv",
    )


def test_anim_gzip(
    tmp_path: str, input_genomes_tiny: Path, gzipped_tiny_example: Path
) -> None:
    """Check ANIm run (gzipped)."""
    tmp_dir = Path(tmp_path) / "ANIm's gzip test 📦"
    tmp_dir.mkdir()
    tmp_db = tmp_dir / "anim's inputs are gzipped.sqlite"
    public_cli.cli_anim(
        database=tmp_db, fasta=gzipped_tiny_example, name="Test Run", create_db=True
    )
    public_cli.export_run(database=tmp_db, outdir=tmp_dir, run_id=1)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "ANIm_identity.tsv",
        tmp_dir / "ANIm_identity.tsv",
    )

    logger = setup_logger(None)
    # Now do it again but with the decompressed files - should reuse calculations:
    session = db_orm.connect_to_db(logger, tmp_db)
    count = session.query(db_orm.Comparison).count()
    session.close()

    public_cli.cli_anim(
        database=tmp_db, fasta=input_genomes_tiny, name="Test Run", create_db=False
    )

    session = db_orm.connect_to_db(logger, tmp_db)
    assert count == session.query(db_orm.Comparison).count()
    session.close()


def test_dnadiff(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
    evil_example: Path,
) -> None:
    """Check dnadiff run (default settings)."""
    caplog.set_level(logging.INFO)
    tmp_dir = Path(tmp_path) / "dnadiff's test 📋"
    tmp_dir.mkdir()
    tmp_db = tmp_dir / "dnadiff test.sqlite"
    # Leaving out name, so can check the default worked
    public_cli.cli_dnadiff(database=tmp_db, fasta=evil_example, create_db=True)
    output = caplog.text
    assert "Database already has 0 of 3²=9 dnadiff comparisons, 9 needed\n" in output
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = session.query(db_orm.Run).one()
    assert run.name == "3 genomes using dnadiff"
    session.close()

    # Now do it again - it should reuse the calculations:
    caplog.clear()
    public_cli.cli_dnadiff(
        database=tmp_db,
        fasta=input_genomes_tiny,
        name="Simple names",
        create_db=False,
    )
    output = caplog.text
    assert "Database already has all 3²=9 dnadiff comparisons\n" in output

    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    # Fuzzy, 0.9963 from dnadiff tool != 0.9962661747 from our code
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "dnadiff_identity.tsv",
        tmp_dir / "dnadiff_identity.tsv",
        atol=5e-5,
    )


def test_dnadiff_gzip(
    tmp_path: str, input_genomes_tiny: Path, gzipped_tiny_example: Path
) -> None:
    """Check dnadiff run (gzipped)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "dnadiff's inputs are gzipped.sqlite"
    # Leaving out name, so can check the default worked
    public_cli.cli_dnadiff(database=tmp_db, fasta=gzipped_tiny_example, create_db=True)
    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    # Fuzzy, 0.9963 from dnadiff tool != 0.9962661747 from our code
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "dnadiff_identity.tsv",
        tmp_dir / "dnadiff_identity.tsv",
        atol=5e-5,
    )
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = session.query(db_orm.Run).one()
    assert run.name == "3 genomes using dnadiff"
    session.close()


def test_anib(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
    evil_example: Path,
) -> None:
    """Check ANIb run (spaces, emoji, etc in filenames)."""
    caplog.set_level(logging.INFO)
    tmp_dir = Path(tmp_path) / "ANIb-test-🎱"  # no spaces! makeblastdb -out breaks
    tmp_dir.mkdir()
    tmp_db = tmp_dir / "anib test.sqlite"
    public_cli.cli_anib(
        database=tmp_db,
        fasta=evil_example,
        name="Spaces etc",
        create_db=True,
        temp=tmp_dir,
    )
    output = caplog.text
    assert "Database already has 0 of 3²=9 ANIb comparisons, 9 needed\n" in output

    # Run it again, nothing to recompute but easier to check output
    caplog.clear()
    public_cli.cli_anib(
        database=tmp_db,
        fasta=input_genomes_tiny,
        name="Simple names",
        create_db=True,
        temp=tmp_dir,
    )
    output = caplog.text
    assert "Database already has all 3²=9 ANIb comparisons\n" in output

    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "ANIb_identity.tsv",
        tmp_dir / "ANIb_identity.tsv",
    )


def test_anib_gzip(
    tmp_path: str, input_genomes_tiny: Path, gzipped_tiny_example: Path
) -> None:
    """Check ANIb run (gzipped)."""
    tmp_dir = Path(tmp_path) / "a:b'c;d-😞"  # no spaces for makeblastdb -out
    tmp_dir.mkdir()
    tmp_db = tmp_dir / "anib's inputs are gzipped.sqlite"
    public_cli.cli_anib(
        database=tmp_db,
        fasta=gzipped_tiny_example,
        name="Test Run",
        create_db=True,
        temp=tmp_dir,
    )

    # The fragment files should match those expected without gzip compression!
    for file in (input_genomes_tiny / "intermediates/ANIb").glob("*.f*"):
        assert filecmp.cmp(file, tmp_dir / file), f"Wrong fragmented FASTA {file.name}"

    # The intermediate TSV files should match too
    for file in (input_genomes_tiny / "intermediates/ANIb").glob("*_vs_*.tsv"):
        assert filecmp.cmp(file, tmp_dir / file), f"Wrong blastn output in {file.name}"

    # Since the matrices are labelled by stem, they should match too:
    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "ANIb_identity.tsv",
        tmp_dir / "ANIb_identity.tsv",
    )


def test_fastani(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
    evil_example: Path,
) -> None:
    """Check fastANI run (spaces, emoji, etc in filenames)."""
    caplog.set_level(logging.INFO)
    tmp_dir = Path(tmp_path) / "fastANI's test 🏎️"
    tmp_dir.mkdir()
    tmp_db = tmp_dir / "fastani's test.sqlite"
    public_cli.cli_fastani(
        database=tmp_db,
        fasta=evil_example,
        name="Spaces etc",
        create_db=True,
        temp=tmp_dir,
    )
    output = caplog.text
    assert "Database already has 0 of 3²=9 fastANI comparisons, 9 needed\n" in output

    # Run it again, nothing to recompute but easier to check output
    caplog.clear()
    public_cli.cli_fastani(
        database=tmp_db,
        fasta=input_genomes_tiny,
        name="Simple names",
        create_db=False,
        temp=tmp_dir,
    )
    output = caplog.text
    assert "Database already has all 3²=9 fastANI comparisons\n" in output

    # Confirm output matches
    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "fastANI_identity.tsv",
        tmp_dir / "fastANI_identity.tsv",
    )


def test_fastani_gzip(tmp_path: str, input_gzip_bacteria: Path) -> None:
    """Check fastANI run (gzipped bacteria)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "fastani's inputs are gzipped.sqlite"
    public_cli.cli_fastani(
        database=tmp_db,
        fasta=input_gzip_bacteria,
        name="Test Run",
        create_db=True,
        temp=tmp_dir,
    )

    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    compare_matrix_files(
        input_gzip_bacteria / "matrices" / "fastANI_identity.tsv",
        tmp_dir / "fastANI_identity.tsv",
    )


def test_sourmash_gzip(tmp_path: str, input_gzip_bacteria: Path) -> None:
    """Check sourmash run (gzipped bacteria)."""
    tmp_dir = Path(tmp_path) / "sourmash's gzip test 🚅"
    tmp_dir.mkdir()
    tmp_db = tmp_dir / "sourmash's inputs are gzipped.sqlite"
    public_cli.cli_sourmash(
        database=tmp_db,
        fasta=input_gzip_bacteria,
        name="Test Run",
        create_db=True,
        cache=tmp_dir,
    )
    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    compare_matrix_files(
        input_gzip_bacteria / "matrices" / "sourmash_identity.tsv",
        tmp_dir / "sourmash_identity.tsv",
    )


def test_sourmash(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
    evil_example: Path,
) -> None:
    """Check sourmash run (default settings except scaled=300)."""
    caplog.set_level(logging.INFO)
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "sourmash test.sqlite"

    public_cli.cli_sourmash(
        database=tmp_db,
        fasta=evil_example,
        name="Spaces etc",
        scaled=300,
        create_db=True,
        cache=tmp_dir,
    )
    output = caplog.text
    assert "Database already has 0 of 3²=9 sourmash comparisons, 9 needed\n" in output

    # Run it again, nothing to recompute but easier to check output
    caplog.clear()
    public_cli.cli_sourmash(
        database=tmp_db,
        fasta=input_genomes_tiny,
        name="Simple names",
        scaled=300,
        create_db=False,
        cache=tmp_dir,
    )
    output = caplog.text
    assert "Database already has all 3²=9 sourmash comparisons\n" in output

    # Confirm output matches
    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "sourmash_identity.tsv",
        tmp_dir / "sourmash_identity.tsv",
    )
    plot_out = tmp_dir / "plots 📊"
    plot_out.mkdir()
    # Should be able to plot run 1 with the spaces and emoji, but get warnings:
    # UserWarning: Glyph 129440 (\N{MICROBE}) missing from font(s) DejaVu Sans.
    public_cli.plot_run(database=tmp_db, outdir=plot_out, run_id=2)
    assert sorted(_.name for _ in plot_out.glob("*")) == sorted(
        [
            f"sourmash_{name}_{kind}.{ext}"
            for name in ("identity", "query_cov", "hadamard", "tANI")
            for kind in ("heatmap", "dist")
            for ext in GRAPHICS_FORMATS
            if not (kind == "dist" and ext == "tsv")
        ]
        + [
            f"sourmash_{name}_scatter.{ext}"
            for name in ("query_cov", "tANI")
            for ext in GRAPHICS_FORMATS
        ]
    )
    for f in sorted(_.name for _ in plot_out.glob("*.tsv")):
        if "_scatter." in f:
            # Currently we don't deliberately sort the scatter points, so they
            # are as per the DB, which in turn depends on how they were imported.
            # For sourmash was predictable, but for branchwater it is random.
            with (
                (plot_out / f).open() as new,
                (input_genomes_tiny / "plots" / Path(f)).open() as old,
            ):
                assert sorted(new.readlines()) == sorted(old.readlines())
        else:
            assert filecmp.cmp(
                plot_out / f, input_genomes_tiny / "plots" / Path(f).name
            ), f


def test_fastani_dups(tmp_path: str) -> None:
    """Check fastANI run (duplicate FASTA inputs)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "fastani-dups.sqlite"
    for name in ("alpha", "beta", "gamma"):
        with (tmp_dir / (name + ".fasta")).open("w") as handle:
            handle.write(">genome\nACGTACGT\n")
    with pytest.raises(SystemExit, match="Multiple genomes with same MD5 checksum"):
        public_cli.cli_fastani(
            database=tmp_db, fasta=tmp_dir, name="Test duplicates fail", create_db=True
        )


def test_resume_partial_fastani(
    caplog: pytest.LogCaptureFixture,
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check list-runs and export-run with mock data including a partial fastANI run."""
    caplog.set_level(logging.INFO)
    tmp_db = Path(tmp_path) / "partial resume.sqlite"
    tool = tools.get_fastani()
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
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
        db_orm.db_genome(logger, session, filename, md5, create=True)

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
    assert " Method  ┃ Done ┃ Null ┃ Miss ┃ Total ┃ Status " in output, output
    assert " fastANI │    8 │    0 │    1 │  9=3² │ Partial " in output, output

    caplog.clear()
    public_cli.resume(database=tmp_db)
    output = caplog.text
    assert "Resuming run-id 1\n" in output, output
    assert "Database already has 8 of 3²=9 fastANI comparisons, 1 needed" in output, (
        output
    )

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " fastANI │    9 │    0 │    0 │  9=3² │ Done " in output, output


def test_resume_partial_anib(
    caplog: pytest.LogCaptureFixture,
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check list-runs and export-run with mock data including a partial ANIb run."""
    caplog.set_level(logging.INFO)
    tmp_db = Path(tmp_path) / "resume.sqlite"
    tool = tools.get_blastn()
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
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
        db_orm.db_genome(logger, session, filename, md5, create=True)

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
    assert " Method ┃ Done ┃ Null ┃ Miss ┃ Total ┃ Status " in output, output
    assert " ANIb   │    8 │    0 │    1 │  9=3² │ Partial " in output, output

    caplog.clear()
    public_cli.resume(database=tmp_db)
    output = caplog.text
    assert "Resuming run-id 1\n" in output, output
    assert "Database already has 8 of 3²=9 ANIb comparisons, 1 needed" in output, output

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " ANIb   │    9 │    0 │    0 │  9=3² │ Done " in output, output


def test_resume_partial_anim(
    caplog: pytest.LogCaptureFixture,
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check list-runs and export-run with mock data including a partial ANIm run."""
    caplog.set_level(logging.INFO)
    tmp_db = Path(tmp_path) / "resume.sqlite"
    tool = tools.get_nucmer()
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
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
        db_orm.db_genome(logger, session, filename, md5, create=True)

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
    assert " Method ┃ Done ┃ Null ┃ Miss ┃ Total ┃ Status " in output, output
    assert " ANIm   │    8 │    0 │    1 │  9=3² │ Partial " in output, output

    caplog.clear()
    public_cli.resume(database=tmp_db)
    output = caplog.text
    assert "Resuming run-id 1\n" in output, output
    assert "Database already has 8 of 3²=9 ANIm comparisons, 1 needed" in output, output

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " ANIm   │    9 │    0 │    0 │  9=3² │ Done " in output, output


def test_resume_partial_sourmash(
    caplog: pytest.LogCaptureFixture,
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check list-runs and export-run with mock data including a partial sourmash run."""
    caplog.set_level(logging.INFO)
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "resume sourmash.sqlite"
    tool = tools.get_sourmash()
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    config = db_orm.db_configuration(
        session,
        "sourmash",
        tool.exe_path.stem,
        tool.version,
        kmersize=31,  # must be 31 to match the sig files in fixtures
        extra="scaled=300",
        create=True,
    )

    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(logger, session, filename, md5, create=True)

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
    assert " Method   ┃ Done ┃ Null ┃ Miss ┃ Total ┃ Status " in output, output
    assert " sourmash │    4 │    0 │    5 │  9=3² │ Partial " in output, output

    caplog.clear()
    public_cli.resume(database=tmp_db, cache=tmp_dir)
    output = caplog.text
    assert "Resuming run-id 1\n" in output, output
    assert "Database already has 4 of 3²=9 sourmash comparisons, 5 needed" in output, (
        output
    )

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " sourmash │    9 │    0 │    0 │  9=3² │ Done " in output, output


def test_resume_dir_gone(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check expected failure trying to resume without the input directory."""
    tmp_db = Path(tmp_path) / "resume.sqlite"
    tool = tools.get_fastani()
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
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
        db_orm.db_genome(logger, session, filename, md5, create=True)

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
        match=r"run-id 1 used input folder /mnt/shared/old, but that is not a directory \(now\).",
    ):
        public_cli.resume(database=tmp_db)


def test_resume_unknown(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check expected failure trying to resume an unknown method."""
    tmp_db = Path(tmp_path) / "resume.sqlite"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
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
        db_orm.db_genome(logger, session, filename, md5, create=True)

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
        match="Unknown method guessing for run-id 1 in .*/resume.sqlite",
    ):
        public_cli.resume(database=tmp_db)


def test_resume_complete(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check resume works for all the methods (using completed runs for speed)."""
    caplog.set_level(logging.INFO)
    tmp_db = Path(tmp_path) / "resume.sqlite"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)

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
            db_orm.db_genome(logger, session, filename, md5, create=True)

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
        caplog.clear()
        public_cli.resume(database=tmp_db)
        output = caplog.text
        assert f"Resuming run-id {index + 1}\n" in output, output
        assert f"Database already has all 3²=9 {method} comparisons" in output, output


def test_resume_fasta_gone(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check resume error handling when a FASTA file is missing."""
    caplog.set_level(logging.INFO)
    tmp_dir = Path(tmp_path)
    tmp_indir = tmp_dir / "input"
    tmp_indir.mkdir()
    tmp_db = tmp_dir / "resume.sqlite"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
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
        db_orm.db_genome(logger, session, filename, md5, create=True)

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
            f"run-id 1 used .*/{Path(missing).name} with MD5 {fasta_to_hash[missing]}"
            " but this FASTA file no longer exists"
        ),
    ):
        public_cli.resume(database=tmp_db)

    # Should work with all the FASTA files, even though now in a different directory
    # to that logged in the genome table (as could happen via an older run):
    (tmp_indir / Path(missing).name).symlink_to(Path(missing))
    caplog.clear()
    public_cli.resume(database=tmp_db)
    output = caplog.text
    assert "Database already has all 3²=9 ANIb comparisons" in output, output

    # Should work even with extra FASTA files, real world use case:
    # Using input directory like /mnt/shared/genomes
    # Start a run when there are 100 genomes (but it was not completed)
    # Add some more genomes, say now 150 genomes
    # Resume the run - it should only operate on the first 100 genomes!
    with (tmp_indir / "extra.fasta").open("w") as handle:
        handle.write(">recently-added-genome\nACGTACGTAGT\n")
    caplog.clear()
    public_cli.resume(database=tmp_db)
    output = caplog.text
    assert "Database already has all 3²=9 ANIb comparisons" in output, output


def test_plot_skip_nulls(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check export-run behaviour when have null values."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "plot_null.sqlite"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
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
        db_orm.db_genome(logger, session, filename, md5, create=True)

    # Record all of the possible comparisons, leaving coverage null
    genomes = list(fasta_to_hash.values())
    for query_hash in genomes:
        for subject_hash in genomes:
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
        cmdline="pyani guessing ...",
        fasta_directory=input_genomes_tiny,
        status="Partial",
        name="Test plotting when some data is null",
        fasta_to_hash=fasta_to_hash,
    )
    session.close()

    plot_out = tmp_dir / "plots"
    caplog.clear()
    caplog.set_level(logging.INFO)
    public_cli.plot_run(database=tmp_db, outdir=plot_out)
    output = caplog.text

    assert f"Output directory {plot_out} does not exist, making it." in output, output
    assert "Cannot plot query_cov as all NA\n" in output, output
    assert "Cannot plot hadamard as all NA\n" in output, output
    assert "Cannot plot tANI as all NA\n" in output, output
    assert f"Wrote {2 * len(GRAPHICS_FORMATS)} images to {plot_out}" in output, output
    assert sorted(_.name for _ in plot_out.glob("*_heatmap.*")) == sorted(
        f"guessing_identity_heatmap.{ext}" for ext in GRAPHICS_FORMATS
    )


def test_plot_bad_nulls(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check export-run behaviour when have null values except on diagonal."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "plot_null.sqlite"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
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
        for filename in sorted(input_genomes_tiny.glob("*.fa*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(logger, session, filename, md5, create=True)

    # Record all of the possible comparisons, leaving coverage null
    genomes = list(fasta_to_hash.values())
    for query_hash in genomes:
        for subject_hash in genomes:
            db_orm.db_comparison(
                session,
                config.configuration_id,
                query_hash,
                subject_hash,
                1.0 if query_hash is subject_hash else None,
                12345,
                cov_query=1.0 if query_hash is subject_hash else None,
            )

    db_orm.add_run(
        session,
        config,
        cmdline="pyani guessing ...",
        fasta_directory=input_genomes_tiny,
        status="Partial",
        name="Test plotting when off-diagonal is null",
        fasta_to_hash=fasta_to_hash,
    )
    session.close()

    plot_out = tmp_dir / "plots"
    plot_out.mkdir()
    caplog.clear()
    public_cli.plot_run(database=tmp_db, outdir=plot_out)
    output = caplog.text
    assert (
        "identity matrix contains 2 nulls (out of 2²=4 guessing comparisons)" in output
    ), output
    assert (
        "query_cov matrix contains 2 nulls (out of 2²=4 guessing comparisons)" in output
    ), output
    assert (
        "hadamard matrix contains 2 nulls (out of 2²=4 guessing comparisons)" in output
    ), output
    assert (
        "tANI matrix contains 2 nulls (out of 2²=4 guessing comparisons)" in output
    ), output


def test_classify_failures(tmp_path: str) -> None:
    """Check classify failures."""
    tmp_dir = Path(tmp_path)

    with pytest.raises(SystemExit, match="Database /does/not/exist does not exist"):
        public_cli.cli_classify(database=Path("/does/not/exist"), outdir=tmp_dir)


def test_classify_warnings(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check classify warnings."""
    tmp_dir = Path(tmp_path)

    # Record only one comparison (self-to-self)
    tmp_db = tmp_dir / "classify_complete.sqlite"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    config = db_orm.db_configuration(
        session, "fastANI", "fastani", "1.2.3", create=True
    )

    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(logger, session, filename, md5, create=True)

    genomes = list(fasta_to_hash.values())
    for query_hash in genomes:
        for subject_hash in genomes:
            db_orm.db_comparison(
                session,
                config.configuration_id,
                query_hash,
                subject_hash,
                cov_query=1.0,
                identity=1.0,
                sim_errors=0,
            )

    db_orm.add_run(
        session,
        config,
        cmdline="pyani fastANI ...",
        fasta_directory=input_genomes_tiny,
        status="Complete",
        name="Test classify when all data present",
        fasta_to_hash=dict(list(fasta_to_hash.items())[2:]),
    )

    db_orm.add_run(
        session,
        config,
        cmdline="pyani fastANI ...",
        fasta_directory=input_genomes_tiny,
        status="Complete",
        name="Test classify when all data present",
        fasta_to_hash=fasta_to_hash,
    )

    caplog.clear()
    public_cli.cli_classify(database=tmp_db, outdir=tmp_dir, run_id=1)
    output = caplog.text
    assert (
        "Run 1 has 1 comparison across 1 genome. Reporting single clique." in output
    ), output

    with (tmp_dir / "fastANI_classify.tsv").open() as handle:
        assert (
            handle.readline()
            == "n_nodes\tmax_cov\tmin_identity\tmax_identity\tmembers\n"
        )
    caplog.clear()
    caplog.set_level(logging.INFO)
    public_cli.cli_classify(database=tmp_db, outdir=tmp_dir, run_id=2, cov_min=1.0)
    output = caplog.text
    assert "All genomes are singletons. No plot can be generated." in output, output
    with (tmp_dir / "fastANI_classify.tsv").open() as handle:
        assert (
            handle.readline()
            == "n_nodes\tmax_cov\tmin_identity\tmax_identity\tmembers\n"
        )


def test_classify_normal(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check working example of classify."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "classify_complete.sqlite"
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    config = db_orm.db_configuration(
        session, "fastANI", "fastani", "1.2.3", create=True
    )

    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(logger, session, filename, md5, create=True)

    # Record all of the possible comparisons
    genomes = list(fasta_to_hash.values())
    cov_increment = 0.88
    for query_hash in genomes:
        for subject_hash in genomes:
            db_orm.db_comparison(
                session,
                config.configuration_id,
                query_hash,
                subject_hash,
                1.0 if query_hash is subject_hash else 0.99,
                12345,
                cov_query=1.0 if query_hash is subject_hash else cov_increment,
            )
            # Increment cov_increment if the hashes are not the same
            if query_hash is not subject_hash:
                cov_increment += 0.01

    db_orm.add_run(
        session,
        config,
        cmdline="pyani fastANI ...",
        fasta_directory=input_genomes_tiny,
        status="Complete",
        name="Test classify when all data present",
        fasta_to_hash=fasta_to_hash,
    )

    caplog.clear()
    caplog.set_level(logging.INFO)
    public_cli.cli_classify(database=tmp_db, outdir=tmp_dir, cov_min=0.9, mode="tANI")
    output = caplog.text
    assert f"Wrote classify output to {tmp_path}" in output, output
    with (tmp_dir / "fastANI_classify.tsv").open() as handle:
        assert handle.readline() == "n_nodes\tmax_cov\tmin_-tANI\tmax_-tANI\tmembers\n"


def test_plot_run_comp(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Plot comparisons with three mock runs."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "⚔️.db"

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    fasta_to_hash = {
        filename: file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(logger, session, filename, md5, create=True)
    genomes = list(fasta_to_hash.values())

    config_a = db_orm.db_configuration(session, "ANIb", "blastn", "0.0", create=True)
    config_b = db_orm.db_configuration(session, "ANIm", "nucmer", "0.0", create=True)
    config_c = db_orm.db_configuration(
        session, "fastANI", "fastani", "0.0", create=True
    )

    for i, config in enumerate((config_a, config_b, config_c)):
        for query_hash in genomes:
            for subject_hash in genomes:
                db_orm.db_comparison(
                    session,
                    config.configuration_id,
                    query_hash,
                    subject_hash,
                    1.0 if query_hash is subject_hash else 0.99 ** (i + 1),
                    12345,
                )
        db_orm.add_run(
            session,
            config,
            cmdline="pyani-plus ...",
            fasta_directory=input_genomes_tiny,
            status="Done",
            name="Trial " + "ABC"[i],
            fasta_to_hash=fasta_to_hash,
        )
    session.commit()
    session.close()

    for cols in (0, 1):
        plot_out = tmp_dir / f"📊{cols}col"
        caplog.clear()
        caplog.set_level(logging.INFO)
        public_cli.plot_run_comp(
            database=tmp_db,
            outdir=plot_out,
            run_ids="1,2,3",
            columns=cols,
        )
        output = caplog.text
        images = len(GRAPHICS_FORMATS)
        if "tsv" in GRAPHICS_FORMATS:
            images -= 1
        assert (
            f"Wrote {images * 2} images to {plot_out}/ANIb_identity_1_vs_*.*\n"
            in output
        ), output

        assert sorted(_.name for _ in plot_out.glob("*")) == sorted(
            [
                f"ANIb_identity_1_{mode}_vs_others.{ext}"
                for ext in GRAPHICS_FORMATS
                for mode in ("scatter", "diff")
                if ext != "tsv"
            ]
            + [f"ANIb_identity_1_vs_{other}.tsv" for other in ("2", "3")]
        )
