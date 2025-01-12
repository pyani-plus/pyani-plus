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
import shutil
from pathlib import Path

import pandas as pd
import pytest

from pyani_plus import db_orm, public_cli, tools
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

    with pytest.raises(SystemExit, match="ERROR: Database contains no runs."):
        public_cli.resume(database=tmp_db)

    with pytest.raises(
        SystemExit,
        match="ERROR: Database has no run-id 1.",
    ):
        public_cli.resume(database=tmp_db, run_id=1)


def test_partial_run(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check list-runs and export-run with mock data including a partial run."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "list-runs.sqlite"
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
        fasta_to_hash=dict(list(fasta_to_hash.items())[1:]),
    )

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 3 analysis runs in " in output, output
    assert " Method  ┃ Done ┃ Null ┃ Miss ┃ Total ┃ Status " in output, output
    assert " fastANI │    0 │    0 │    0 │  0=0² │ Empty " in output, output
    assert " fastANI │    4 │    0 │    5 │  9=3² │ Partial " in output, output
    assert " fastANI │    4 │    0 │    0 │  4=2² │ Done " in output, output

    # Unlike a typical method calculation, we have not triggered
    # .cache_comparisons() yet, so that will happen in export_run.
    public_cli.export_run(database=tmp_db, run_id=3, outdir=tmp_path, label="md5")
    output = capsys.readouterr().out
    assert f"Wrote matrices to {tmp_path}" in output, output
    with (tmp_dir / "fastANI_identity.tsv").open() as handle:
        assert (
            handle.readline()
            == "\t5584c7029328dc48d33f95f0a78f7e57\t78975d5144a1cd12e98898d573cf6536\n"
        )

    public_cli.export_run(database=tmp_db, run_id=3, outdir=tmp_path, label="stem")
    output = capsys.readouterr().out
    assert f"Wrote matrices to {tmp_path}" in output, output
    with (tmp_dir / "fastANI_identity.tsv").open() as handle:
        assert handle.readline() == "\tMGV-GENOME-0266457\tOP073605\n"

    public_cli.export_run(database=tmp_db, run_id=3, outdir=tmp_path, label="filename")
    output = capsys.readouterr().out
    assert f"Wrote matrices to {tmp_path}" in output, output
    with (tmp_dir / "fastANI_identity.tsv").open() as handle:
        assert handle.readline() == "\tMGV-GENOME-0266457.fna\tOP073605.fasta\n"

    # By construction run 2 is partial, only 4 of 9 matrix entries are
    # defined - this should fail
    with pytest.raises(
        SystemExit,
        match=("ERROR: run-id 2 has only 4 of 3²=9 comparisons, 5 needed"),
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


def test_export_run_failures(tmp_path: str) -> None:
    """Check export run failures."""
    tmp_dir = Path(tmp_path)

    with pytest.raises(
        SystemExit, match="ERROR: Database /does/not/exist does not exist"
    ):
        public_cli.export_run(database=Path("/does/not/exist"), outdir=tmp_dir)

    with pytest.raises(
        SystemExit, match="ERROR: Output directory /does/not/exist does not exist"
    ):
        public_cli.export_run(database=":memory:", outdir=Path("/does/not/exist/"))

    tmp_db = tmp_dir / "export.sqlite"
    session = db_orm.connect_to_db(tmp_db)
    with pytest.raises(SystemExit, match="ERROR: Database contains no runs."):
        public_cli.export_run(database=tmp_db, outdir=tmp_dir)

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
    with pytest.raises(SystemExit, match="ERROR: Database has no run-id 3."):
        public_cli.export_run(database=tmp_db, outdir=tmp_dir, run_id=3)
    with pytest.raises(
        SystemExit,
        match="ERROR: run-id 1 has no comparisons",
    ):
        public_cli.export_run(database=tmp_db, outdir=tmp_dir, run_id=1)
    # Should default to latest run, run-id 2
    with pytest.raises(
        SystemExit,
        match="ERROR: run-id 2 has no comparisons",
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
    session = db_orm.connect_to_db(tmp_db)
    config = db_orm.db_configuration(
        session, "fastANI", "fastani", "1.2.3", create=True
    )
    fasta_to_hash = {fasta: file_md5sum(fasta) for fasta in tmp_fasta.glob("*.fa*")}
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
            )
    session.commit()
    session.close()

    with pytest.raises(
        SystemExit,
        match="ERROR: Duplicate filename stems, consider using MD5 labelling.",
    ):
        public_cli.export_run(database=tmp_db, outdir=tmp_dir)

    with pytest.raises(
        SystemExit,
        match="ERROR: Duplicate filename stems, consider using MD5 labelling.",
    ):
        public_cli.plot_run(database=tmp_db, outdir=tmp_dir)


def test_plot_run_failures(tmp_path: str) -> None:
    """Check plot run failures."""
    tmp_dir = Path(tmp_path)

    with pytest.raises(
        SystemExit, match="ERROR: Database /does/not/exist does not exist"
    ):
        public_cli.plot_run(database=Path("/does/not/exist"), outdir=tmp_dir)

    with pytest.raises(
        SystemExit, match="ERROR: Output directory /does/not/exist does not exist"
    ):
        public_cli.plot_run(database=":memory:", outdir=Path("/does/not/exist/"))

    tmp_db = tmp_dir / "export.sqlite"
    session = db_orm.connect_to_db(tmp_db)
    with pytest.raises(SystemExit, match="ERROR: Database contains no runs."):
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
    with pytest.raises(SystemExit, match="ERROR: Database has no run-id 3."):
        public_cli.plot_run(database=tmp_db, outdir=tmp_dir, run_id=3)
    with pytest.raises(
        SystemExit,
        match="ERROR: run-id 1 has no comparisons",
    ):
        public_cli.plot_run(database=tmp_db, outdir=tmp_dir, run_id=1)
    # Should default to latest run, run-id 2
    with pytest.raises(
        SystemExit,
        match="ERROR: run-id 2 has no comparisons",
    ):
        public_cli.plot_run(database=tmp_db, outdir=tmp_dir)


def test_anim(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check ANIm run."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "example.sqlite"
    public_cli.anim(
        database=tmp_db, fasta=input_genomes_tiny, name="Test Run", create_db=True
    )
    public_cli.export_run(database=tmp_db, outdir=tmp_dir, run_id=1)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "ANIm_identity.tsv",
        tmp_dir / "ANIm_identity.tsv",
    )

    # Now do it again - it should reuse the calculations:
    session = db_orm.connect_to_db(tmp_db)
    count = session.query(db_orm.Comparison).count()
    session.close()

    public_cli.anim(
        database=tmp_db, fasta=input_genomes_tiny, name="Test Run", create_db=False
    )

    session = db_orm.connect_to_db(tmp_db)
    assert count == session.query(db_orm.Comparison).count()
    session.close()


def test_anim_gzip(
    tmp_path: str, input_genomes_tiny: Path, gzipped_tiny_example: Path
) -> None:
    """Check ANIm run (gzipped)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "example.sqlite"
    public_cli.anim(
        database=tmp_db, fasta=gzipped_tiny_example, name="Test Run", create_db=True
    )
    public_cli.export_run(database=tmp_db, outdir=tmp_dir, run_id=1)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "ANIm_identity.tsv",
        tmp_dir / "ANIm_identity.tsv",
    )

    # Now do it again but with the decompressed files - it should reuse the calculations:
    session = db_orm.connect_to_db(tmp_db)
    count = session.query(db_orm.Comparison).count()
    session.close()

    public_cli.anim(
        database=tmp_db, fasta=input_genomes_tiny, name="Test Run", create_db=False
    )
    session = db_orm.connect_to_db(tmp_db)
    assert count == session.query(db_orm.Comparison).count()
    session.close()


def test_dnadiff(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check dnadiff run (default settings)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "example.sqlite"
    # Leaving out name, so can check the default worked
    public_cli.dnadiff(database=tmp_db, fasta=input_genomes_tiny, create_db=True)
    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    # Fuzzy, 0.9963 from dnadiff tool != 0.9962661747 from our code
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "dnadiff_identity.tsv",
        tmp_dir / "dnadiff_identity.tsv",
        atol=5e-5,
    )
    session = db_orm.connect_to_db(tmp_db)
    run = session.query(db_orm.Run).one()
    assert run.name == "3 genomes using dnadiff"
    session.close()


def test_dnadiff_gzip(
    tmp_path: str, input_genomes_tiny: Path, gzipped_tiny_example: Path
) -> None:
    """Check dnadiff run (gzipped)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "example.sqlite"
    # Leaving out name, so can check the default worked
    public_cli.dnadiff(database=tmp_db, fasta=gzipped_tiny_example, create_db=True)
    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    # Fuzzy, 0.9963 from dnadiff tool != 0.9962661747 from our code
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "dnadiff_identity.tsv",
        tmp_dir / "dnadiff_identity.tsv",
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


def test_anib_gzip(
    tmp_path: str, input_genomes_tiny: Path, gzipped_tiny_example: Path
) -> None:
    """Check ANIb run (gzipped)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "example.sqlite"
    public_cli.anib(
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


def test_fastani_gzip(tmp_path: str, input_gzip_bacteria: Path) -> None:
    """Check fastANI run (gzipped bacteria)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "example.sqlite"
    public_cli.fastani(
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
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "example.sqlite"
    public_cli.sourmash(
        database=tmp_db,
        fasta=input_gzip_bacteria,
        name="Test Run",
        create_db=True,
    )
    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    compare_matrix_files(
        input_gzip_bacteria / "matrices" / "sourmash_identity.tsv",
        tmp_dir / "sourmash_identity.tsv",
    )


def test_sourmash(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check sourmash run (default settings except scaled=300)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "example.sqlite"
    public_cli.sourmash(
        database=tmp_db,
        fasta=input_genomes_tiny,
        name="Test Run",
        scaled=300,
        create_db=True,
    )
    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "sourmash_identity.tsv",
        tmp_dir / "sourmash_identity.tsv",
    )
    plot_out = tmp_dir / "plots"
    plot_out.mkdir()
    public_cli.plot_run(database=tmp_db, outdir=plot_out)
    assert sorted(_.name for _ in plot_out.glob("*")) == [
        "sourmash_hadamard.jpg",
        "sourmash_hadamard.pdf",
        "sourmash_hadamard.png",
        "sourmash_hadamard.svg",
        "sourmash_hadamard.tsv",
        "sourmash_identity.jpg",
        "sourmash_identity.pdf",
        "sourmash_identity.png",
        "sourmash_identity.svg",
        "sourmash_identity.tsv",
        "sourmash_query_cov.jpg",
        "sourmash_query_cov.pdf",
        "sourmash_query_cov.png",
        "sourmash_query_cov.svg",
        "sourmash_query_cov.tsv",
    ]


def test_branchwater(tmp_path: str, input_genomes_tiny: Path) -> None:
    """Check sourmash run (default settings except scaled=300)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "example.sqlite"
    public_cli.branchwater(
        database=tmp_db,
        fasta=input_genomes_tiny,
        name="Test Run",
        scaled=300,
        create_db=True,
    )
    public_cli.export_run(database=tmp_db, outdir=tmp_dir)
    # Should match the sourmash output (but computed quicker)
    compare_matrix_files(
        input_genomes_tiny / "matrices" / "sourmash_identity.tsv",
        tmp_dir / "branchwater_identity.tsv",
    )


def test_fastani_dups(tmp_path: str) -> None:
    """Check fastANI run (duplicate FASTA inputs)."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "example.sqlite"
    for name in ("alpha", "beta", "gamma"):
        with (tmp_dir / (name + ".fasta")).open("w") as handle:
            handle.write(">genome\nACGTACGT\n")
    with pytest.raises(
        SystemExit, match="ERROR - Multiple genomes with same MD5 checksum"
    ):
        public_cli.fastani(
            database=tmp_db, fasta=tmp_dir, name="Test duplicates fail", create_db=True
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
    assert " Method  ┃ Done ┃ Null ┃ Miss ┃ Total ┃ Status " in output, output
    assert " fastANI │    8 │    0 │    1 │  9=3² │ Partial " in output, output

    public_cli.resume(database=tmp_db)
    output = capsys.readouterr().out
    assert "Resuming run-id 1\n" in output, output
    assert "Database already has 8 of 3²=9 fastANI comparisons, 1 needed" in output, (
        output
    )

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " fastANI │    9 │    0 │    0 │  9=3² │ Done " in output, output


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
    assert " Method ┃ Done ┃ Null ┃ Miss ┃ Total ┃ Status " in output, output
    assert " ANIb   │    8 │    0 │    1 │  9=3² │ Partial " in output, output

    public_cli.resume(database=tmp_db)
    output = capsys.readouterr().out
    assert "Resuming run-id 1\n" in output, output
    assert "Database already has 8 of 3²=9 ANIb comparisons, 1 needed" in output, output

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " ANIb   │    9 │    0 │    0 │  9=3² │ Done " in output, output


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
    assert " Method ┃ Done ┃ Null ┃ Miss ┃ Total ┃ Status " in output, output
    assert " ANIm   │    8 │    0 │    1 │  9=3² │ Partial " in output, output

    public_cli.resume(database=tmp_db)
    output = capsys.readouterr().out
    assert "Resuming run-id 1\n" in output, output
    assert "Database already has 8 of 3²=9 ANIm comparisons, 1 needed" in output, output

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " ANIm   │    9 │    0 │    0 │  9=3² │ Done " in output, output


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
    assert " Method   ┃ Done ┃ Null ┃ Miss ┃ Total ┃ Status " in output, output
    assert " sourmash │    4 │    0 │    5 │  9=3² │ Partial " in output, output

    public_cli.resume(database=tmp_db)
    output = capsys.readouterr().out
    assert "Resuming run-id 1\n" in output, output
    assert "Database already has 4 of 3²=9 sourmash comparisons, 5 needed" in output, (
        output
    )

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    assert " sourmash │    9 │    0 │    0 │  9=3² │ Done " in output, output


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
    assert " Method   ┃ Done ┃ Null ┃ Miss ┃ Total ┃ Status " in output, output
    assert " branchw… │    4 │    0 │    5 │  9=3² │ Partial " in output, output

    public_cli.resume(database=tmp_db)
    output = capsys.readouterr().out
    assert "Resuming run-id 1\n" in output, output
    assert (
        "Database already has 4 of 3²=9 branchwater comparisons, 5 needed" in output
    ), output

    public_cli.list_runs(database=tmp_db)
    output = capsys.readouterr().out
    assert " 1 analysis runs in " in output, output
    # Note the layout shifted from the above as the status column
    # can be narrower giving a character more for the method:
    assert " Method    ┃ Done ┃ Null ┃ Miss ┃ Total ┃ Status" in output, output
    assert " branchwa… │    9 │    0 │    0 │  9=3² │ Done " in output, output


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
        assert f"Resuming run-id {index + 1}\n" in output, output
        assert f"Database already has all 3²=9 {method} comparisons" in output, output


def test_resume_fasta_gone(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check resume error handling when a FASTA file is missing."""
    tmp_dir = Path(tmp_path)
    tmp_indir = tmp_dir / "input"
    tmp_indir.mkdir()
    tmp_db = tmp_dir / "resume.sqlite"
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
    assert "Database already has all 3²=9 ANIb comparisons" in output, output

    # Should work even with extra FASTA files, real world use case:
    # Using input directory like /mnt/shared/genomes
    # Start a run when there are 100 genomes (but it was not completed)
    # Add some more genomes, say now 150 genomes
    # Resume the run - it should only operate on the first 100 genomes!
    with (tmp_indir / "extra.fasta").open("w") as handle:
        handle.write(">recently-added-genome\nACGTACGTAGT\n")
    public_cli.resume(database=tmp_db)
    output = capsys.readouterr().out
    assert "Database already has all 3²=9 ANIb comparisons" in output, output


def test_plot_skip_nulls(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check export-run behaviour when have null values."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "plot_null.sqlite"
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
    plot_out.mkdir()
    public_cli.plot_run(database=tmp_db, outdir=plot_out)

    stdout, stderr = capsys.readouterr()
    assert (
        "WARNING: Cannot plot query_cov as matrix contains 9 nulls (out of 3²=9 guessing comparisons)\n"
        in stderr
    ), stderr
    assert (
        "WARNING: Cannot plot hadamard as matrix contains 9 nulls (out of 3²=9 guessing comparisons)\n"
        in stderr
    ), stderr
    assert "Wrote 1 heatmaps" in stdout
    assert sorted(_.name for _ in plot_out.glob("*")) == [
        "guessing_identity.jpg",
        "guessing_identity.pdf",
        "guessing_identity.png",
        "guessing_identity.svg",
        "guessing_identity.tsv",
    ]


def test_plot_bad_nulls(
    capsys: pytest.CaptureFixture[str], tmp_path: str, input_genomes_tiny: Path
) -> None:
    """Check export-run behaviour when have null values except on diagonal."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "plot_null.sqlite"
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
        for filename in sorted(input_genomes_tiny.glob("*.fa*"))
    }
    for filename, md5 in fasta_to_hash.items():
        db_orm.db_genome(session, filename, md5, create=True)

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

    with pytest.raises(
        SystemExit,
        match=(r"ERROR: Unable to plot any heatmaps \(check for nulls\)"),
    ):
        public_cli.plot_run(database=tmp_db, outdir=plot_out)

    stderr = capsys.readouterr().err
    assert (
        "WARNING: Cannot plot identity as matrix contains 2 nulls (out of 2²=4 guessing comparisons)\n"
        in stderr
    ), stderr
    assert (
        "WARNING: Cannot plot query_cov as matrix contains 2 nulls (out of 2²=4 guessing comparisons)\n"
        in stderr
    ), stderr
    assert (
        "WARNING: Cannot plot hadamard as matrix contains 2 nulls (out of 2²=4 guessing comparisons)\n"
        in stderr
    ), stderr
