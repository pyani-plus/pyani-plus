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
"""Tests for interrupting pyANI-plus jobs.

These tests are intended to be run from the repository root using:

pytest -v
"""

import contextlib
import signal
import subprocess
import time
from pathlib import Path

import pytest
from Bio.SeqIO.FastaIO import SimpleFastaParser

from pyani_plus import db_orm, private_cli, tools

GENOMES = 100


@pytest.fixture(scope="session")
def large_set_of_bacterial_chunks(
    tmp_path_factory: pytest.TempPathFactory, input_bacteria: Path
) -> Path:
    """Make a directory of FASTA files, each a chunk of bacteria."""
    fasta_dir = tmp_path_factory.mktemp(f"{GENOMES}_faked_bacteria")
    # Make input dataset of fragments of a bacteria
    with (input_bacteria / "NC_010338.fna").open() as handle:
        title, seq = next(SimpleFastaParser(handle))
    for i in range(GENOMES):
        offset = i * 1000
        with (fasta_dir / f"genome{i}.fasta").open("w") as handle:
            handle.write(f">genome{i}\n{seq[offset:offset+500000]}\n")
    return fasta_dir


# Note using capfd rather than capsys due to subprocesses
def test_compute_column_sigint_anib(
    capfd: pytest.CaptureFixture[str],
    tmp_path: str,
    large_set_of_bacterial_chunks: Path,
) -> None:
    """Check compute_column with SIGINT and ANIb."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "sigint.sqlite"
    assert not tmp_db.is_file()

    tmp_tool = tmp_dir / "working"
    tmp_tool.mkdir()

    tmp_input = tmp_dir / "input"
    tmp_input.symlink_to(large_set_of_bacterial_chunks)

    tool = tools.get_blastn()
    with contextlib.chdir(tmp_dir):
        private_cli.log_run(
            fasta=Path("input/"),  # Not absolute but relative to test that works
            database=tmp_db,
            cmdline="pyani-plus anib ...",
            status="Testing",
            name="Testing SIGINT",
            method="ANIb",
            program=tool.exe_path.stem,
            version=tool.version,
            fragsize=1000,
            create_db=True,
        )
    output = capfd.readouterr().out
    assert output.endswith("Run identifier 1\n")

    with subprocess.Popen(
        [
            ".pyani-plus-private-cli",
            "compute-column",
            "--quiet",
            "--database",
            str(tmp_db),
            "--run-id",
            "1",
            "--subject",
            # Compute the final column which would not do its partial logging
            # until the very end (unless interrupted):
            str(GENOMES - 1),
            "--temp",
            "working/",  # Not absolute but relative to test that works
        ],
        cwd=tmp_dir,
    ) as process:
        time.sleep(15)
        process.send_signal(signal.SIGINT)
        process.wait()
        rc = process.returncode
    output = capfd.readouterr().err

    session = db_orm.connect_to_db(tmp_db)
    done = session.query(db_orm.Comparison).count()
    genomes = session.query(db_orm.Genome).count()

    if done < genomes:
        assert "Interrupted, will attempt to log" in output, output
    else:
        assert (
            0 < done < genomes
        ), f"Expected partial ANIb progress, not {done}/{genomes} and return {rc}"
    assert (
        rc == 0
    ), f"Expecting return 0 on clean interruption, got {rc} with {done} done"


def test_compute_column_sigint_dnadiff(
    capfd: pytest.CaptureFixture[str],
    tmp_path: str,
    large_set_of_bacterial_chunks: Path,
) -> None:
    """Check compute_column with SIGINT and dnadiff."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "sigint.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_nucmer()
    private_cli.log_run(
        fasta=large_set_of_bacterial_chunks,
        database=tmp_db,
        cmdline="pyani-plus dnadiff ...",
        status="Testing",
        name="Testing SIGINT",
        method="dnadiff",
        program=tool.exe_path.stem,
        version=tool.version,
        create_db=True,
    )
    output = capfd.readouterr().out
    assert output.endswith("Run identifier 1\n")

    with subprocess.Popen(
        [
            ".pyani-plus-private-cli",
            "compute-column",
            "--quiet",
            "--database",
            str(tmp_db),
            "--run-id",
            "1",
            "--subject",
            # Compute the final column which would not do its partial logging
            # until the very end (unless interrupted):
            str(GENOMES - 1),
        ],
        cwd=tmp_dir,
    ) as process:
        time.sleep(10)
        process.send_signal(signal.SIGINT)
        process.wait()
        rc = process.returncode
    output = capfd.readouterr().err

    session = db_orm.connect_to_db(tmp_db)
    done = session.query(db_orm.Comparison).count()
    genomes = session.query(db_orm.Genome).count()

    if done < genomes:
        assert "Interrupted, will attempt to log" in output, output
    else:
        assert (
            0 < done < genomes
        ), f"Expected partial dnadiff progress, not {done}/{genomes} and return {rc}"
    assert (
        rc == 0
    ), f"Expecting return 0 on clean interruption, got {rc} with {done} done"


def test_compute_column_sigint_anim(
    capfd: pytest.CaptureFixture[str],
    tmp_path: str,
    large_set_of_bacterial_chunks: Path,
) -> None:
    """Check compute_column with SIGINT and ANIm."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "sigint.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_nucmer()
    private_cli.log_run(
        fasta=large_set_of_bacterial_chunks,
        database=tmp_db,
        cmdline="pyani-plus anim ...",
        status="Testing",
        name="Testing SIGINT",
        method="ANIm",
        mode="mum",
        program=tool.exe_path.stem,
        version=tool.version,
        create_db=True,
    )
    output = capfd.readouterr().out
    assert output.endswith("Run identifier 1\n")

    with subprocess.Popen(
        [
            ".pyani-plus-private-cli",
            "compute-column",
            "--quiet",
            "--database",
            str(tmp_db),
            "--run-id",
            "1",
            "--subject",
            # Compute the final column which would not do its partial logging
            # until the very end (unless interrupted):
            str(GENOMES - 1),
        ],
        cwd=tmp_dir,
    ) as process:
        time.sleep(8)
        process.send_signal(signal.SIGINT)
        process.wait()
        rc = process.returncode
    output = capfd.readouterr().err

    session = db_orm.connect_to_db(tmp_db)
    done = session.query(db_orm.Comparison).count()
    genomes = session.query(db_orm.Genome).count()

    if done < genomes:
        assert "Interrupted, will attempt to log" in output, output
    else:
        assert (
            0 < done < genomes
        ), f"Expected partial ANIb progress, not {done}/{genomes} and return {rc}"
    assert (
        rc == 0
    ), f"Expecting return 0 on clean interruption, got {rc} with {done} done"