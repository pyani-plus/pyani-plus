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
"""Tests for the sourmash implementation.

These tests are intended to be run from the repository root using:

make test
"""

import subprocess
from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, tools, utils
from pyani_plus.methods import method_sourmash

from . import get_matrix_entry


def test_compare_parsing(sourmash_targets_compare_indir: Path) -> None:
    """Check parsing of sourmash compare .csv file."""
    assert (
        method_sourmash.parse_compare(
            sourmash_targets_compare_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv"
        )
        == 0.997900938305757  # noqa: PLR2004
    )


def test_missing_db(
    tmp_path: str, input_genomes_tiny: Path, sourmash_targets_compare_indir: Path
) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.log_sourmash(
            database=tmp_db,
            # These are for the comparison table
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            compare=sourmash_targets_compare_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv",
        )


def test_bad_query_or_subject(
    tmp_path: str, input_genomes_tiny: Path, sourmash_targets_compare_indir: Path
) -> None:
    """Mismatch between query or subject FASTA in sourmash compare output and commandline."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --query-fasta .*/MGV-GENOME-0266457.fna"
            " but query in sourmash compare filename was MGV-GENOME-0264574"
        ),
    ):
        private_cli.log_sourmash(
            database=tmp_db,
            query_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
            compare=sourmash_targets_compare_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv",
        )

    with pytest.raises(
        SystemExit,
        match=(
            "ERROR: Given --subject-fasta .*/MGV-GENOME-0264574.fas"
            " but subject in sourmash compare filename was MGV-GENOME-0266457"
        ),
    ):
        private_cli.log_sourmash(
            database=tmp_db,
            query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            subject_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
            compare=sourmash_targets_compare_indir
            / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv",
        )


def test_logging_sourmash(
    tmp_path: str,
    input_genomes_tiny: Path,
    sourmash_targets_compare_indir: Path,
    dir_sourmash_results: Path,
) -> None:
    """Check can log a sourmash comparison to DB."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    tool = tools.get_sourmash()

    private_cli.log_configuration(
        database=tmp_db,
        method="sourmash",
        program=tool.exe_path.stem,
        version=tool.version,
        mode=method_sourmash.MODE,
        kmersize=method_sourmash.KMER_SIZE,
        extra="scaled=" + str(method_sourmash.SCALED),
        create_db=True,
    )
    private_cli.log_genome(
        database=tmp_db,
        fasta=[
            input_genomes_tiny / "MGV-GENOME-0264574.fas",
            input_genomes_tiny / "MGV-GENOME-0266457.fna",
        ],
    )

    private_cli.log_sourmash(
        database=tmp_db,
        query_fasta=input_genomes_tiny / "MGV-GENOME-0264574.fas",
        subject_fasta=input_genomes_tiny / "MGV-GENOME-0266457.fna",
        compare=sourmash_targets_compare_indir
        / "MGV-GENOME-0264574_vs_MGV-GENOME-0266457.csv",
    )

    # Check the recorded comparison values
    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 1
    comp = session.query(db_orm.Comparison).one()
    pytest.approx(
        comp.identity,
        get_matrix_entry(
            dir_sourmash_results / "matrix_identity.tsv",
            comp.query_hash,
            comp.subject_hash,
        ),
    )
    session.close()
    tmp_db.unlink()


def test_prepare_compute_sourmash(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Confirm prepare and compute with sourmash."""
    tmp_db = Path(tmp_path) / "done.db"
    cache = Path(tmp_path) / "cache"
    tool = tools.get_sourmash()
    session = db_orm.connect_to_db(tmp_db)
    config = db_orm.db_configuration(
        session,
        "sourmash",
        tool.exe_path.stem,
        tool.version,
        create=True,
    )

    fasta_to_hash = {
        filename: utils.file_md5sum(filename)
        for filename in sorted(input_genomes_tiny.glob("*.f*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(session, filename, md5, create=True)

    db_orm.add_run(
        session,
        config,
        cmdline="pyani-plus sourmash ...",
        fasta_directory=input_genomes_tiny,
        status="Pending",
        name="An initially empty mock run",
        fasta_to_hash=fasta_to_hash,
    )
    session.close()

    # First with no cache,
    assert not cache.is_dir()
    with pytest.raises(
        SystemExit,
        match=f"ERROR: Specified cache directory {cache} does not exist",
    ):
        private_cli.prepare(database=tmp_db, run_id=1, cache=cache)

    with pytest.raises(
        SystemExit,
        match="ERROR: sourmash needs prepared files but cache /.*/cache directory does not exist",
    ):
        private_cli.compute(database=tmp_db, run_id=1, cache=cache)

    # Now with empty cache,
    cache.mkdir()

    # Now run prepare to add the sig file to the cache,
    # but the run configuration is broken
    with pytest.raises(
        SystemExit,
        match="ERROR: sourmash requires a k-mer size, default is 31",
    ):
        private_cli.prepare(database=tmp_db, run_id=1, cache=cache)

    # Let's fix the k-mer config
    config = session.query(db_orm.Configuration).one()
    config.kmersize = 75
    session.commit()

    with pytest.raises(
        SystemExit,
        match="ERROR: sourmash requires scaled or num, default is scaled=1000",
    ):
        private_cli.prepare(database=tmp_db, run_id=1, cache=cache)

    # Let's fix the extra config
    config = session.query(db_orm.Configuration).one()
    config.extra = "scaled=250"  # must be small for virus, test fixures use 300
    session.commit()

    with pytest.raises(
        SystemExit,
        match="ERROR: sourmash requires mode, default is max-containment",
    ):
        private_cli.compute(database=tmp_db, run_id=1, cache=cache)

    # Let's fix the mode config
    config = session.query(db_orm.Configuration).one()
    config.mode = "containment"
    session.commit()

    with pytest.raises(
        subprocess.CalledProcessError,
        # This should fail due to missing .sig file
        # ValueError: Error while reading signature
        match="Command '.*' returned non-zero exit status 1.",
    ):
        private_cli.compute(database=tmp_db, run_id=1, cache=cache)

    # Now the prepare step should work (and make all three .sig files),
    private_cli.prepare(database=tmp_db, run_id=1, cache=cache)
    assert len(list(cache.glob("*.sig"))) == 3  # noqa: PLR2004

    # And now this should work too, note will accept task=0 or N to mean same:
    private_cli.compute(database=tmp_db, run_id=1, cache=cache, task=3, parts=3)
    assert session.query(db_orm.Comparison).count() == 3  # noqa: PLR2004
    private_cli.compute(database=tmp_db, run_id=1, cache=cache, task=0, parts=3)
    assert session.query(db_orm.Comparison).count() == 3  # noqa: PLR2004

    private_cli.compute(database=tmp_db, run_id=1, cache=cache, task=1, parts=3)
    assert session.query(db_orm.Comparison).count() == 6  # noqa: PLR2004
    private_cli.compute(database=tmp_db, run_id=1, cache=cache, task=2, parts=3)
    assert session.query(db_orm.Comparison).count() == 9  # noqa: PLR2004
