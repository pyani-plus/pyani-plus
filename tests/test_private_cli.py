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
from multiprocessing import Pool
from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli


def test_log_configuration(tmp_path: str) -> None:
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


def test_log_run(tmp_path: str) -> None:
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
            fasta=[],  # list
            # Config
            method="guessing",
            program="guestimate",
            version="0.1.2beta3",
            fragsize=100,
            kmersize=51,
            # Misc
            create_db=False,
        )

    # This time create it
    private_cli.log_run(
        database=tmp_db,
        # Run
        cmdline="pyani_plus run ...",
        name="Guess Run",
        status="Completed",
        fasta=[],  # list
        # Config
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        # Misc
        create_db=True,
    )

    tmp_db.unlink()


def test_log_comparison_no_db(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Confirm can create a mock DB using log-comparison alone."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist, but not using --create-db"):
        private_cli.log_comparison(
            database=tmp_db,
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            identity=0.96,
            aln_length=12345,
            method="guessing",
            program="guestimate",
            version="0.1.2beta3",
            fragsize=100,
            kmersize=51,
            sim_errors=1,
            cov_query=0.98,
            cov_subject=0.98,
            create_db=False,
        )


def test_log_comparison_duplicate(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Confirm no error logging comparison twice."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    private_cli.log_comparison(
        database=tmp_db,
        query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
        subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
        identity=0.96,
        aln_length=12345,
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        sim_errors=1,
        cov_query=0.98,
        cov_subject=0.98,
        create_db=True,
    )

    private_cli.log_comparison(
        database=tmp_db,
        query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
        subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
        identity=0.955,  # different!
        aln_length=12345,
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
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


def test_log_comparison_serial(tmp_path: str, input_genomes_tiny: Path) -> None:
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

    fasta = list(input_genomes_tiny.glob("*.fna"))  # subset of folder (2)
    private_cli.log_genome(
        database=tmp_db,
        fasta=fasta,
        create_db=False,
    )

    # Could at this point log the run with status=started (or similar),
    # but will need a mechanism to return the run ID and use it to update
    # the table row at the end...

    for query in fasta:
        for subject in fasta:
            private_cli.log_comparison(
                database=tmp_db,
                query_fasta=query,
                subject_fasta=subject,
                identity=1.0 if query == subject else 0.96,
                aln_length=12345,
                method="guessing",
                program="guestimate",
                version="0.1.2beta3",
                fragsize=100,
                kmersize=51,
                sim_errors=1,
                cov_query=0.98,
                cov_subject=0.98,
                create_db=False,
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
        fasta=fasta,  # list
        # Config
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        # Misc
        create_db=False,
    )

    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == len(fasta) ** 2


def test_log_comparison_parallel(tmp_path: str, input_genomes_tiny: Path) -> None:
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

    fasta = list(input_genomes_tiny.glob("*.f*"))
    pool = Pool(3)
    for filename in fasta:
        # Deliberately add each file multiple times to try to clash
        for _ in range(3):
            pool.apply_async(
                private_cli.log_genome,
                [],
                {
                    "database": tmp_db,
                    "fasta": [filename],
                    "create_db": False,
                },
            )
    pool.close()
    pool.join()

    private_cli.log_genome(
        database=tmp_db,
        fasta=fasta,
        create_db=False,
    )

    # Could at this point log the run with status=started (or similar),
    # but will need a mechanism to return the run ID and use it to update
    # the table row at the end...

    assert tmp_db.is_file()
    tasks = [
        {
            "database": tmp_db,
            "query_fasta": query,
            "subject_fasta": subject,
            "identity": 1.0 if query == subject else 0.96,
            "aln_length": 12345,
            "method": "guessing",
            "program": "guestimate",
            "version": "0.1.2beta3",
            "fragsize": 100,
            "kmersize": 51,
            "sim_errors": 1,
            "cov_query": 0.98,
            "cov_subject": 0.98,
            "create_db": False,
        }
        for query in fasta
        for subject in fasta
    ]

    pool = Pool(len(fasta) ** 2)
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
        fasta=fasta,  # list
        # Config
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        # Misc
        create_db=False,
    )

    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == len(fasta) ** 2


def test_fragment_fasta(
    tmp_path: str, input_genomes_tiny: Path, anib_fragments: Path
) -> None:
    """Confirm fragmenting FASTA files (for ANIb) works."""
    fasta = input_genomes_tiny.glob("*.f*")
    out_dir = Path(tmp_path)
    private_cli.fragment_fasta(fasta, out_dir)

    old_frags = [anib_fragments / (f.stem + "-fragments.fna") for f in fasta]
    new_frags = [out_dir / (f.stem + "-fragments.fna") for f in fasta]
    for old_file, new_file in zip(old_frags, new_frags, strict=True):
        assert filecmp.cmp(old_file, new_file), f"Wrong output in {new_file}"


def test_fragment_fasta_bad_args(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check error handling for fragmenting FASTA files (for ANIb)."""
    fasta = input_genomes_tiny.glob("*.f*")
    out_dir = Path(tmp_path)

    with pytest.raises(
        SystemExit, match="ERROR: outdir /does/not/exist should be a directory"
    ):
        private_cli.fragment_fasta(fasta, outdir=Path("/does/not/exist"))

    with pytest.raises(
        ValueError, match="Cannot fragment /does/not/exist.fasta, file not found"
    ):
        private_cli.fragment_fasta([Path("/does/not/exist.fasta")], out_dir)

    with pytest.raises(
        TypeError,
        match="Expected a Path object in list of FASTA files, got 'some/example.fasta'",
    ):
        private_cli.fragment_fasta(["some/example.fasta"], out_dir)
