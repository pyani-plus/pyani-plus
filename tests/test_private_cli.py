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
"""Tests for the pyani_plus/private_cli.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

import multiprocessing
import sys
from pathlib import Path

import pytest
from sqlalchemy.exc import NoResultFound

from pyani_plus import db_orm, private_cli, tools


def test_log_configuration(capsys: pytest.CaptureFixture[str], tmp_path: str) -> None:
    """Confirm can create a new empty database via log-configuration."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist, but not using --create-db"):
        private_cli.log_configuration(
            tmp_db,
            method="guessing",
            program="guestimate",
            version="0.1.2beta3",
            fragsize=100,
            kmersize=51,
            create_db=False,
        )

    # This time create it
    private_cli.log_configuration(
        tmp_db,
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Configuration identifier 1\n")

    # This time should already be a DB there
    private_cli.log_configuration(
        tmp_db,
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=75,
        kmersize=31,
        create_db=False,
    )
    output = capsys.readouterr().out
    assert output.endswith("Configuration identifier 2\n")

    tmp_db.unlink()


def test_log_genome(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Confirm can create a new empty database via log-genome."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist, but not using --create-db"):
        private_cli.log_genome(
            database=tmp_db,
            fasta=list(
                input_genomes_tiny.glob("*.fasta")  # subset of folder
            ),
        )

    # This time create it
    private_cli.log_genome(
        database=tmp_db,
        fasta=list(
            input_genomes_tiny.glob("*.fasta")  # subset of folder
        ),
        create_db=True,
    )


def test_log_run(capsys: pytest.CaptureFixture[str], tmp_path: str) -> None:
    """Confirm can create a new empty DB via log-run."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist, but not using --create-db"):
        private_cli.log_run(
            database=tmp_db,
            # Run
            cmdline="pyani_plus run ...",
            name="Guess Run",
            status="Completed",
            fasta=Path("/does/not/exist/"),
            # Config
            method="guessing",
            program="guestimate",
            version="0.1.2beta3",
            fragsize=100,
            kmersize=51,
            # Misc
            create_db=False,
        )

    with pytest.raises(SystemExit, match="ERROR: No FASTA input genomes under"):
        private_cli.log_run(
            database=tmp_db,
            # Run
            cmdline="pyani_plus run ...",
            name="Guess Run",
            status="Completed",
            fasta=tmp_path,
            # Config
            method="guessing",
            program="guestimate",
            version="0.1.2beta3",
            fragsize=100,
            kmersize=51,
            # Misc
            create_db=True,
        )

    with (Path(tmp_path) / "example.fasta").open("w") as handle:
        handle.write(">Tiny\nACGTACGTTA\n")

    # This time create it
    private_cli.log_run(
        database=tmp_db,
        # Run
        cmdline="pyani_plus run ...",
        name="Guess Run",
        status="Completed",
        fasta=tmp_path,
        # Config
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        # Misc
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    tmp_db.unlink()


def test_log_comparison_no_db(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Confirm log-comparison fails if DB is missing."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.log_comparison(
            database=tmp_db,
            config_id=1,
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            identity=0.96,
            aln_length=12345,
            sim_errors=1,
            cov_query=0.98,
            cov_subject=0.98,
        )


def test_log_comparison_no_config(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Confirm log-comparison fails if config is missing."""
    tmp_db = Path(tmp_path) / "empty.sqlite"
    assert not tmp_db.is_file()
    session = db_orm.connect_to_db(tmp_db)
    session.commit()
    session.close()

    with pytest.raises(
        SystemExit, match="empty.sqlite does not contain configuration_id=1"
    ):
        private_cli.log_comparison(
            database=tmp_db,
            config_id=1,
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            identity=0.96,
            aln_length=12345,
            sim_errors=1,
            cov_query=0.98,
            cov_subject=0.98,
        )


def test_log_comparison_duplicate(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Confirm no error logging comparison twice."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    private_cli.log_configuration(
        database=tmp_db,
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Configuration identifier 1\n")

    private_cli.log_genome(
        database=tmp_db,
        fasta=[
            input_genomes_tiny / "MGV-GENOME-0264574.fas",
            input_genomes_tiny / "MGV-GENOME-0266457.fna",
        ],
    )

    private_cli.log_comparison(
        database=tmp_db,
        config_id=1,
        query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
        subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
        identity=0.96,
        aln_length=12345,
        sim_errors=1,
        cov_query=0.98,
        cov_subject=0.98,
    )

    private_cli.log_comparison(
        database=tmp_db,
        config_id=1,
        query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
        subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
        identity=0.955,  # different!
        aln_length=12345,
        sim_errors=1,
        cov_query=0.98,
        cov_subject=0.98,
    )

    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 1
    comp = session.query(db_orm.Comparison).one()
    # first value should not be replaced:
    assert comp.identity == 0.96  # noqa: PLR2004
    session.close()
    tmp_db.unlink()


def test_log_comparison_serial_and_skip_process_genomes(
    caplog: pytest.LogCaptureFixture,
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Confirm can create a mock DB using log-comparison etc. sequentially.

    Also test process-genomes detects a complete run and aborts.
    """
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "serial.sqlite"
    assert not tmp_db.is_file()

    private_cli.log_configuration(
        tmp_db,
        method="sourmash",
        program="guestimate",
        version="0.1.2beta3",
        kmersize=51,
        extra="scaled=1234",
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Configuration identifier 1\n")

    fasta = list(input_genomes_tiny.glob("*.f*"))
    private_cli.log_genome(
        database=tmp_db,
        fasta=fasta,
    )

    # Could at this point log the run with status=started (or similar),
    # but will need a mechanism to return the run ID and use it to update
    # the table row at the end...

    for query in fasta:
        for subject in fasta:
            private_cli.log_comparison(
                database=tmp_db,
                config_id=1,
                query_fasta=query,
                subject_fasta=subject,
                identity=1.0 if query == subject else 0.96,
                aln_length=12345,
                sim_errors=1,
                cov_query=0.98,
                cov_subject=0.98,
            )

    # Can now log the run with status=completed
    # Or, if we already logged it with status=started, would need to update
    # the existing run table entry with the cached matrices and completed status
    private_cli.log_run(
        database=tmp_db,
        # Run
        cmdline="pyani_plus run ...",
        name="Guess Run",
        status="Completed",
        fasta=input_genomes_tiny,
        # Config
        method="sourmash",
        program="guestimate",
        version="0.1.2beta3",
        kmersize=51,
        extra="scaled=1234",
        create_db=False,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == len(fasta) ** 2
    assert session.query(db_orm.Configuration).count() == 1

    caplog.clear()
    private_cli.prepare_genomes(database=tmp_db, run_id=1, cache=tmp_dir)
    output = caplog.text
    assert "Skipping preparation, run already has all 9=3² pairwise values" in output


def test_log_comparison_parallel(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Confirm can create a mock DB using log-comparison etc. in parallel."""
    tmp_db = Path(tmp_path) / "parallel.sqlite"
    assert not tmp_db.is_file()

    private_cli.log_configuration(
        tmp_db,
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Configuration identifier 1\n")

    fasta = list(input_genomes_tiny.glob("*.f*"))
    # Avoid implicit fork, should match the defaults on Python 3.14 onwards:
    pool = multiprocessing.get_context(  # type:ignore [attr-defined]
        "spawn" if sys.platform == "darwin" else "forkserver"
    ).Pool(3)
    for filename in fasta:
        # Deliberately add each file multiple times to try to clash
        for _ in range(3):
            pool.apply_async(
                private_cli.log_genome,
                [],
                {
                    "database": tmp_db,
                    "fasta": [filename],
                },
            )
    pool.close()
    pool.join()

    # Could at this point log the run with status=started (or similar),
    # but will need a mechanism to return the run ID and use it to update
    # the table row at the end...

    assert tmp_db.is_file()
    tasks = [
        {
            "database": tmp_db,
            "config_id": 1,
            "query_fasta": query,
            "subject_fasta": subject,
            "identity": 1.0 if query == subject else 0.96,
            "aln_length": 12345,
            "sim_errors": 1,
            "cov_query": 0.98,
            "cov_subject": 0.98,
        }
        for query in fasta
        for subject in fasta
    ]

    # Avoid implicit fork, should match the defaults on Python 3.14 onwards:
    pool = multiprocessing.get_context(  # type:ignore [attr-defined]
        "spawn" if sys.platform == "darwin" else "forkserver"
    ).Pool(len(fasta) ** 2)
    for kwargs in tasks:
        pool.apply_async(private_cli.log_comparison, [], kwargs)
    pool.close()
    pool.join()

    # Can now log the run with status=completed
    # Or, if we already logged it with status=started, would need to update
    # the existing run table entry with the cached matrices and completed status
    private_cli.log_run(
        database=tmp_db,
        # Run
        cmdline="pyani_plus run ...",
        name="Guess Run",
        status="Completed",
        fasta=input_genomes_tiny,
        # Config
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        # Misc
        create_db=False,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == len(fasta) ** 2


def test_validate_cache(tmp_path: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """Check expected error handling in validate_cache."""
    tmp_dir = Path(tmp_path)
    monkeypatch.chdir(tmp_dir)

    assert not Path(".cache").is_dir()

    with pytest.raises(
        SystemExit,
        match="ERROR: Specified cache directory /does/not/exist does not exist",
    ):
        private_cli.validate_cache(Path("/does/not/exist"))

    default = private_cli.validate_cache(None, create_default=False)
    assert default == Path(".cache")
    assert not default.is_dir()
    with pytest.raises(
        SystemExit, match="Default cache directory .cache does not exist."
    ):
        private_cli.validate_cache(None, create_default=False, require=True)

    default = private_cli.validate_cache(None, create_default=True)
    assert default == Path(".cache")
    assert default.is_dir()


def test_missing_db(tmp_path: str) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.prepare_genomes(
            database=tmp_db,
            run_id=1,
        )

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.compute_column(
            database=tmp_db,
            run_id=1,
            subject=1,
        )


def test_prepare_genomes_bad_args(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check error handling in prepare-genomes."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad.sqlite"
    assert not tmp_db.is_file()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus sourmash ...",
        status="Testing",
        name="Testing compute-column",
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        create_db=True,
    )

    with pytest.raises(
        SystemExit,
        match="ERROR: Unknown method guessing, check tool version?",
    ):
        private_cli.prepare_genomes(database=tmp_db, run_id=1, cache=tmp_dir)


def test_compute_column_bad_args(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check how compute_column handles bad run ID or subject."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus guessing ...",
        status="Testing",
        name="Testing compute-column",
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    # If this was the public API, should handle it more gracefully:
    with pytest.raises(
        NoResultFound,
        match="No row was found when one was required",
    ):
        private_cli.compute_column(
            database=tmp_db,
            run_id=2,
            subject="XXXX",
        )

    with pytest.raises(
        SystemExit,
        match="ERROR: Unknown method guessing for run-id 1 in .*/new.sqlite",
    ):
        private_cli.compute_column(
            database=tmp_db,
            run_id=1,
            subject="1",
        )

    with pytest.raises(
        SystemExit,
        match="ERROR: Did not recognise 'XXXX' as an MD5 hash, filename, or column number in run-id 1",
    ):
        private_cli.compute_column(
            database=tmp_db,
            run_id=1,
            subject="XXXX",
        )

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Single column should be in range 1 to 3,"
            " or for some methods 0 meaning all columns, but not -1"
        ),
    ):
        private_cli.compute_column(
            database=tmp_db,
            run_id=1,
            subject="-1",
        )

    with pytest.raises(
        SystemExit,
        match="ERROR: All columns currently only implemented for sourmash",
    ):
        private_cli.compute_column(
            database=tmp_db,
            run_id=1,
            subject="0",
        )


def test_compute_column_bad_anib(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check how compute_column handles bad ANIb settings."""
    tmp_db = Path(tmp_path) / "anib.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_blastn()
    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus guessing ...",
        status="Testing",
        name="Testing compute-column",
        method="ANIb",
        program=tool.exe_path.stem,
        version=tool.version,
        # fragsize=...,  <-- missing!
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    with pytest.raises(
        SystemExit,
        match="ERROR: ANIb run-id 1 is missing fragsize parameter",
    ):
        private_cli.compute_column(
            database=tmp_db,
            run_id=1,
            subject="1",
        )


def test_compute_column_bad_anim(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check how compute_column handles bad ANIm settings."""
    tmp_db = Path(tmp_path) / "anim.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_nucmer()
    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus guessing ...",
        status="Testing",
        name="Testing compute-column",
        method="ANIm",
        program=tool.exe_path.stem,
        version=tool.version,
        # mode=...,  <-- missing!
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    with pytest.raises(
        SystemExit,
        match="ERROR: ANIm run-id 1 is missing mode parameter",
    ):
        private_cli.compute_column(
            database=tmp_db,
            run_id=1,
            subject="1",
        )


def test_compute_column_bad_fastani(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check how compute_column handles bad fastani settings."""
    tmp_db = Path(tmp_path) / "fastani.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_fastani()
    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus guessing ...",
        status="Testing",
        name="Testing compute-column",
        method="fastANI",
        program=tool.exe_path.stem,
        version=tool.version,
        # fragsize=...,  <-- missing!
        create_db=True,
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 1\n")

    with pytest.raises(
        SystemExit,
        match="ERROR: fastANI run-id 1 is missing fragsize parameter",
    ):
        private_cli.compute_column(
            database=tmp_db,
            run_id=1,
            subject="1",
        )

    tool = tools.get_fastani()
    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus guessing ...",
        status="Testing",
        name="Testing compute-column",
        method="fastANI",
        program=tool.exe_path.stem,
        version=tool.version,
        fragsize=1000,
        # kmersize=...,  <-- missing!
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 2\n")

    with pytest.raises(
        SystemExit,
        match="ERROR: fastANI run-id 2 is missing kmersize parameter",
    ):
        private_cli.compute_column(
            database=tmp_db,
            run_id=2,
            subject="1",
        )

    tool = tools.get_fastani()
    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus guessing ...",
        status="Testing",
        name="Testing compute-column",
        method="fastANI",
        program=tool.exe_path.stem,
        version=tool.version,
        fragsize=1000,
        kmersize=9,
        # minmatch=..., <-- missing!
    )
    output = capsys.readouterr().out
    assert output.endswith("Run identifier 3\n")

    with pytest.raises(
        SystemExit,
        match="ERROR: fastANI run-id 3 is missing minmatch parameter",
    ):
        private_cli.compute_column(
            database=tmp_db,
            run_id=3,
            subject="1",
        )


def test_compute_column_fastani(
    caplog: pytest.LogCaptureFixture,
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check compute_column with valid args (using fastANI)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "new.sqlite"
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

    private_cli.compute_column(
        database=tmp_db,
        run_id=1,
        subject="1",  # here passing column number
        temp=tmp_dir,
    )
    output = caplog.text
    assert "Calling fastANI for 3 queries vs 5584c7029328dc48d33f95f0a78f7e57" in output

    # This time should skip any computation:
    caplog.clear()
    private_cli.compute_column(
        database=tmp_db,
        run_id=1,
        subject="5584c7029328dc48d33f95f0a78f7e57",  # here passing hash
    )
    output = caplog.text
    assert (
        "No fastANI comparisons needed against 5584c7029328dc48d33f95f0a78f7e57"
        in output
    )

    # Again, should skip any computation:
    caplog.clear()
    private_cli.compute_column(
        database=tmp_db,
        run_id=1,
        subject="OP073605.fasta",  # here passing filename
    )
    output = caplog.text
    assert (
        "No fastANI comparisons needed against 5584c7029328dc48d33f95f0a78f7e57"
        in output
    )

    # Don't need prepare-genomes with fastANI, but should get this message
    caplog.clear()
    private_cli.prepare_genomes(database=tmp_db, run_id=1, quiet=False)
    output = caplog.text
    assert "No per-genome preparation required for fastANI" in output
