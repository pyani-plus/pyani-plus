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


def test_log_comparison_serial(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Confirm can create a mock DB using log-comparison etc. sequentially."""
    tmp_db = Path(tmp_path) / "serial.sqlite"
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

    fasta = list(input_genomes_tiny.glob("*.fna"))  # subset of folder
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


def test_log_wrong_config(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Confirm can't run/log a comparison to another method."""
    tmp_db = Path(tmp_path) / "guestimate.sqlite"
    assert not tmp_db.is_file()
    faked = Path(tmp_path) / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.tmp"
    with faked.open("w") as handle:
        handle.write("# Faked\n")

    private_cli.log_run(
        database=tmp_db,
        # Run
        cmdline="pyani_plus run ...",
        name="Guess Run",
        status="Empty",
        fasta=input_genomes_tiny,
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

    with pytest.raises(SystemExit, match="ERROR: Run-id 1 expected guessing results"):
        private_cli.log_anim(
            tmp_db,
            run_id=1,
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            deltafilter=faked,
        )

    with pytest.raises(SystemExit, match="ERROR: Run-id 1 expected guessing results"):
        private_cli.log_dnadiff(
            tmp_db,
            run_id=1,
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            mcoords=faked,
            qdiff=faked,
        )

    with pytest.raises(SystemExit, match="ERROR: Run-id 1 expected guessing results"):
        private_cli.log_sourmash(
            tmp_db,
            run_id=1,
            compare=faked,
        )

    with pytest.raises(SystemExit, match="ERROR: Run-id 1 expected guessing results"):
        private_cli.log_branchwater(
            tmp_db,
            run_id=1,
            manysearch=faked,
        )


def test_compute_column_missing_db(tmp_path: str) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.compute_column(
            database=tmp_db,
            run_id=1,
            subject=0,
        )


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

    with pytest.raises(  # noqa: PT012
        SystemExit,
        match="ERROR: Unknown method guessing for run-id 1 in .*/new.sqlite",
    ):
        private_cli.compute_column(
            database=tmp_db,
            run_id=1,
            subject="3",
        )
        # This has to be within the context manager to work:
        output = capsys.readouterr().out
        assert "INFO: Treating subject N as 0 (first column)" in output

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
        match="ERROR: Column should be in range 0 to 3, not -1",
    ):
        private_cli.compute_column(
            database=tmp_db,
            run_id=1,
            subject="-1",
        )

    with pytest.raises(
        SystemExit,
        match="ERROR: Column should be in range 0 to 3, not 4",
    ):
        private_cli.compute_column(
            database=tmp_db,
            run_id=1,
            subject="4",
        )


def test_compute_column_fastani(
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
        subject="0",  # here passing column number
        temp=tmp_dir,
    )
    output = capsys.readouterr().out
    assert (
        "Calling fastANI for 3 queries vs 5584c7029328dc48d33f95f0a78f7e57\n" in output
    )

    # This time should skip any computation:
    private_cli.compute_column(
        database=tmp_db,
        run_id=1,
        subject="5584c7029328dc48d33f95f0a78f7e57",  # here passing hash
    )
    output = capsys.readouterr().out
    assert (
        "INFO: No fastANI comparisons needed against 5584c7029328dc48d33f95f0a78f7e57\n"
        in output
    )

    # Again, should skip any computation:
    private_cli.compute_column(
        database=tmp_db,
        run_id=1,
        subject="OP073605.fasta",  # here passing filename
    )
    output = capsys.readouterr().out
    assert (
        "INFO: No fastANI comparisons needed against 5584c7029328dc48d33f95f0a78f7e57\n"
        in output
    )
