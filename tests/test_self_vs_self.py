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
"""Tests for self-vs-self comparisons which do not get 100% identity scores.

These tests are intended to be run from the repository root using:

pytest -v
"""

import logging
import tempfile
from collections.abc import Callable
from pathlib import Path

import pytest

from pyani_plus import db_orm, public_cli, setup_logger


def do_self_compare(
    caplog: pytest.LogCaptureFixture,
    fasta_dir: Path,
    method: Callable,
    **kwargs: Path,
) -> db_orm.Comparison:
    """Run a self-vs-self ANI method and return the single comparison."""
    assert len(list(fasta_dir.glob("*.f*"))) == 1
    with tempfile.NamedTemporaryFile(suffix=".db") as tmp_db:
        caplog.set_level(logging.INFO)
        caplog.clear()
        method(
            database=tmp_db.name,
            fasta=fasta_dir,
            name="Self-vs-self",
            create_db=True,
            **kwargs,
        )
        output = caplog.text
        assert " run setup with 1 genomes in database\n" in output

        logger = setup_logger(None)
        session = db_orm.connect_to_db(logger, tmp_db.name)
        comp = session.query(db_orm.Comparison).one()
        session.close()
        return comp


def test_self_vs_self_anim(caplog: pytest.LogCaptureFixture, tmp_path: str) -> None:
    """Check a self comparison where ANIm only 99.6307% identity.

    This is a real contig from GCA_001829655.1_ASM182965v1_genomic.fna
    whole genome shotgun sequence of Sulfurimonas sp.

    It has one run of 28N, but nothing else odd by eye.
    """
    tmp_dir = Path(tmp_path)
    seq_dir = tmp_dir / "fasta"
    seq_dir.mkdir()
    (seq_dir / "MIBY01000005.fasta").symlink_to(
        Path("tests/fixtures/MIBY01000005.fasta").resolve()
    )
    assert len(list(seq_dir.glob("*.f*"))) == 1

    comp = do_self_compare(caplog, seq_dir, public_cli.cli_anim)
    assert comp.identity == 0.9963070429965708, comp  # noqa: PLR2004

    comp = do_self_compare(caplog, seq_dir, public_cli.cli_dnadiff)
    assert comp.identity == 1.0, comp

    comp = do_self_compare(caplog, seq_dir, public_cli.cli_anib)
    assert comp.identity == 1.0, comp

    comp = do_self_compare(caplog, seq_dir, public_cli.cli_fastani)
    assert comp.identity == 1.0, comp

    comp = do_self_compare(caplog, seq_dir, public_cli.cli_sourmash, cache=tmp_dir)
    assert comp.identity == 1.0, comp


def test_self_vs_self_fastani(caplog: pytest.LogCaptureFixture, tmp_path: str) -> None:
    """Check a self comparison where fastANI only 99.9953% identity.

    This is a real contig from GCA_001829655.1_ASM182965v1_genomic.fna
    whole genome shotgun sequence of Sulfurimonas sp.
    """
    tmp_dir = Path(tmp_path)
    seq_dir = tmp_dir / "fasta"
    seq_dir.mkdir()
    (seq_dir / "MIBY01000011.fasta").symlink_to(
        Path("tests/fixtures/MIBY01000011.fasta").resolve()
    )
    assert len(list(seq_dir.glob("*.f*"))) == 1

    comp = do_self_compare(caplog, seq_dir, public_cli.cli_anim)
    assert comp.identity == 1.0, comp

    comp = do_self_compare(caplog, seq_dir, public_cli.cli_dnadiff)
    assert comp.identity == 1.0, comp

    comp = do_self_compare(caplog, seq_dir, public_cli.cli_anib)
    assert comp.identity == 1.0, comp

    comp = do_self_compare(caplog, seq_dir, public_cli.cli_fastani)
    assert comp.identity == 0.999953, comp  # noqa: PLR2004

    comp = do_self_compare(caplog, seq_dir, public_cli.cli_sourmash, cache=tmp_dir)
    assert comp.identity == 1.0, comp
