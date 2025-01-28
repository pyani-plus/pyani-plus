# The MIT License
#
# Copyright (c) 2025 University of Strathclyde
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
"""Tests for importing external alignments and computing their ANI.

These tests are intended to be run from the repository root using:

make test
"""

# Required to support pytest automated testing
from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, public_cli
from pyani_plus.utils import file_md5sum


def test_missing_db(tmp_path: str) -> None:
    """Check expected error when DB does not exist."""
    tmp_db = Path(tmp_path) / "new.sqlite"
    assert not tmp_db.is_file()

    with pytest.raises(SystemExit, match="does not exist"):
        public_cli.external_alignment(
            database=tmp_db,
            fasta=Path("/dev/null"),  # won't get as far as opening this
            alignment=Path("/dev/null"),
        )

    with pytest.raises(SystemExit, match="does not exist"):
        private_cli.log_external_alignment(
            database=tmp_db,
            run_id=1,
        )


def test_simple_mock_alignment_stem(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Mock alignments using filename stem naming.

    These are 9 or 10bp, giving 10bp long pairwise alignments.
    """
    tmp_db = Path(tmp_path) / "stems.db"
    assert not tmp_db.is_file()

    tmp_alignment = Path(tmp_path) / "stems.fasta"

    with tmp_alignment.open("w") as handle:
        # Listing the fragments in MD5 order to match the matrix in DB
        handle.write("""\
>OP073605 mock 10bp fragment for 5584c7029328dc48d33f95f0a78f7e57
AACC-GGTTTT
>MGV-GENOME-0264574 mock 9bp fragment for 689d3fd6881db36b5e08329cf23cecdd
AACC-GG-TTT
>MGV-GENOME-0266457 mock 10bp fragment for 78975d5144a1cd12e98898d573cf6536
AACC-GGATTT
    """)

    public_cli.external_alignment(
        input_genomes_tiny, tmp_db, create_db=True, alignment=tmp_alignment
    )
    output = capsys.readouterr().out
    assert "\nexternal-alignment run setup with 3 genomes" in output, output

    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 9  # noqa: PLR2004
    # Consider the diagonals - the aligned fragments are 10, 9, and 10 bp
    # and their genomes are 57793, 39253, and 39594 bp. Therefore the
    # self-vs-self coverage values are 10/57793, 9/39253, and 10/39594.
    #
    # Consider this pair (1st and 2nd entries):
    #
    # 5584c7029328dc48d33f95f0a78f7e57 AACC-GGTTTT query length 10, genome len 57793
    # 689d3fd6881db36b5e08329cf23cecdd AACC-GG-TTT subject length 9, genome len 39253
    #
    # Equivalent to this:
    #
    # 5584c7029328dc48d33f95f0a78f7e57 AACCGGTTTT query length 10
    # 689d3fd6881db36b5e08329cf23cecdd AACCGG-TTT subject length 9
    #
    # 9 matches 0 mismatch in alignment of gapped-length 10, so identity 9/10 = 0.9
    # However query coverage (9+0)/57793, and subject cover (9+0)/39253
    #
    # So identity [1.0 0.9 ...]    query coverage [10/57793 9/57793 ...     ]
    #             [0.9 1.0 ...]                   [9/39253  9/39253 ...     ]
    #             [... ... 1.0]                   [...      ...     10/39594]
    #
    # Now consider this pair (1st and 3rd entries):
    #
    # 5584c7029328dc48d33f95f0a78f7e57 AACC-GGTTTT query length 10, genome len 57793
    # 78975d5144a1cd12e98898d573cf6536 AACC-GGATTT subject length 10, genome len 39594
    #
    # Equivalent to this:
    #
    # 5584c7029328dc48d33f95f0a78f7e57 AACCGGTTTT
    # 78975d5144a1cd12e98898d573cf6536 AACCGGATTT
    #
    # 9 matches 1 mismatch in an alignment of length 10, identity 9/10=0.9
    # Query coverage (9+1)/57793, subject coverage (9+1)/39594
    #
    # So identity [1.0 ... 0.9]    query coverage [10/57793  ...    10/57793]
    #             [... 1.0 ...]                   [  ...    9/39253   ...   ]
    #             [0.9 ... 1.0]                   [10/39594  ...    10/39594]
    #
    # Finally, consider this pair (2nd and 3rd entries):
    #
    # 689d3fd6881db36b5e08329cf23cecdd AACC-GG-TTT query length 9, genome len 39253
    # 78975d5144a1cd12e98898d573cf6536 AACC-GGATTT subject length 10, genome len 39594
    #
    # Equivalent to this:
    #
    # 689d3fd6881db36b5e08329cf23cecdd AACCGG-TTT query length 9
    # 78975d5144a1cd12e98898d573cf6536 AACCGGATTT subject length 10
    #
    # 9 matches 0 mismatch in alignment of length 9, so identity 9/10 = 0.9
    # However query coverage (9+0)/57793, subject cover (9+0)/39594
    #
    # So identity [1.0 ... ...]    query coverage [10/57793   ...     ...   ]
    #             [... 1.0 0.9]                   [  ...    9/39253  9/39253]
    #             [... 0.9 1.0]                   [  ...    9/39594 10/39594]
    #
    # Overall:
    #
    # So identity [1.0 0.9 0.9]    query coverage [10/57793 9/57793 10/57793]
    #             [0.9 1.0 0.9]                   [ 9/39253 9/39253  9/39253]
    #             [0.9 0.9 1.0]                   [10/39594 9/39594 10/39594]
    #
    run = db_orm.load_run(session, 1, check_complete=True)
    assert run.df_identity == (
        '{"columns":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"index":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"data":[[1.0,0.9,0.9],[0.9,1.0,0.9],[0.9,0.9,1.0]]}'
    )
    assert run.df_cov_query == (
        '{"columns":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"index":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"data":['
        f"[{10 / 57793:.10f},{9 / 57793:.10f},{10 / 57793:.10f}],"
        f"[{9 / 39253:.10f},{9 / 39253:.10f},{9 / 39253:.10f}],"
        f"[{10 / 39594:.10f},{9 / 39594:.10f},{10 / 39594:.10f}"
        "]]}"
    )


def test_simple_mock_alignment_md5(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Mock alignments using MD5 naming.

    These are 4 or 5bp, giving 5bp long pairwise alignments.
    """
    tmp_db = Path(tmp_path) / "filenames.db"
    assert not tmp_db.is_file()

    tmp_alignment = Path(tmp_path) / "filenames.fasta"
    with tmp_alignment.open("w") as handle:
        handle.write("""\
>5584c7029328dc48d33f95f0a78f7e57 aka OP073605.fasta  mock 5bp fragment
CGGTT
>689d3fd6881db36b5e08329cf23cecdd aka MGV-GENOME-0264574.fas mock 4bp fragment
CGG-T
>78975d5144a1cd12e98898d573cf6536 aka MGV-GENOME-0266457.fna mock 5bp fragment
CGGAT
    """)

    public_cli.external_alignment(
        input_genomes_tiny, tmp_db, create_db=True, alignment=tmp_alignment, label="md5"
    )
    output = capsys.readouterr().out
    assert "\nexternal-alignment run setup with 3 genomes" in output, output

    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 9  # noqa: PLR2004
    run = db_orm.load_run(session, 1, check_complete=True)
    assert run.df_identity == (
        '{"columns":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"index":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"data":[[1.0,0.8,0.8],[0.8,1.0,0.8],[0.8,0.8,1.0]]}'
    )


def test_simple_mock_alignment_filename(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Mock alignments using filename naming.

    These are 4 or 5bp, giving 5bp long pairwise alignments.
    """
    tmp_db = Path(tmp_path) / "filenames.db"
    assert not tmp_db.is_file()

    tmp_alignment = Path(tmp_path) / "filenames.fasta"
    with tmp_alignment.open("w") as handle:
        handle.write("""\
>OP073605.fasta mock 5bp fragment
-C-GGTT
>MGV-GENOME-0266457.fna mock 5bp fragment
-C-GGAT
>MGV-GENOME-0264574.fas mock 4bp fragment
-C-GG-T
    """)

    public_cli.external_alignment(
        input_genomes_tiny,
        tmp_db,
        create_db=True,
        alignment=tmp_alignment,
        label="filename",
    )
    output = capsys.readouterr().out
    assert "\nexternal-alignment run setup with 3 genomes" in output, output

    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 9  # noqa: PLR2004
    run = db_orm.load_run(session, 1, check_complete=True)
    assert run.df_identity == (
        '{"columns":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"index":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"data":[[1.0,0.8,0.8],[0.8,1.0,0.8],[0.8,0.8,1.0]]}'
    )


def test_resume(
    capsys: pytest.CaptureFixture[str],
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Resume importing an external alignment."""
    tmp_db = Path(tmp_path) / "stems.db"
    assert not tmp_db.is_file()

    tmp_alignment = Path(tmp_path) / "stems.fasta"
    with tmp_alignment.open("w") as handle:
        handle.write("""\
>OP073605 mock 10bp fragment
AACC-GGTTTT
>MGV-GENOME-0266457 mock 10bp fragment
AACC-GGATTT
>MGV-GENOME-0264574 mock 9bp fragment
AACC-GG-TTT
    """)

    md5 = file_md5sum(tmp_alignment)

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus external-alignment ...",
        status="Pending",
        name="Testing resume",
        method="external-alignment",
        program="",
        version="",
        extra=f"md5={md5};label=stem;alignment={tmp_alignment}",
        create_db=True,
    )

    public_cli.resume(tmp_db)
    output = capsys.readouterr().out
    assert (
        "Database already has 0 of 3²=9 external-alignment comparisons, 9 needed"
        in output
    ), output

    session = db_orm.connect_to_db(tmp_db)
    assert session.query(db_orm.Comparison).count() == 9  # noqa: PLR2004
    run = db_orm.load_run(session, 1, check_complete=True)
    assert run.df_identity == (
        '{"columns":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"index":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"data":[[1.0,0.9,0.9],[0.9,1.0,0.9],[0.9,0.9,1.0]]}'
    )


def test_bad_resume(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Resume with unexpected tool version information."""
    tmp_db = Path(tmp_path) / "stems.db"
    assert not tmp_db.is_file()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus external-alignment ...",
        status="Pending",
        name="Testing resume",
        method="external-alignment",
        program="should-be-blank",
        version="1.0",
        extra="md5=...;label=stem;alignment=...",
        create_db=True,
    )

    with pytest.raises(
        SystemExit,
        match="ERROR: We expect no tool information, but run-id 1 used should-be-blank version 1.0 instead.",
    ):
        public_cli.resume(tmp_db)


def test_bad_mapping(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Mock alignments using bad filename naming.

    These are 4 or 5bp, giving 5bp long pairwise alignments.
    """
    tmp_db = Path(tmp_path) / "filenames.db"
    assert not tmp_db.is_file()

    tmp_alignment = Path(tmp_path) / "filenames.fasta"
    with tmp_alignment.open("w") as handle:
        handle.write("""\
>OP073605.fasta mock 5bp fragment
CGGTT
>MGV-GENOME-0266457.fna mock 5bp fragment
CGGAT
>MGV-GENOME-0264574.fna with wrong filename
CGG-T
    """)

    with pytest.raises(SystemExit, match="Snakemake workflow failed"):
        public_cli.external_alignment(
            input_genomes_tiny,
            tmp_db,
            create_db=True,
            alignment=tmp_alignment,
            label="filename",
        )


def test_wrong_method(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check log-external-alignment rejects a run for another method."""
    tmp_db = Path(tmp_path) / "wrong-method.sqlite"
    assert not tmp_db.is_file()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus guessing ...",
        status="Testing",
        name="Testing compute-column",
        method="guessing",
        program="guestimate",
        version="1.0",
        create_db=True,
    )

    with pytest.raises(
        SystemExit,
        match="ERROR: Run-id 1 expected guessing results",
    ):
        private_cli.log_external_alignment(database=tmp_db, run_id=1)


def test_no_config(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check how log-external-alignment handles bad settings."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad.sqlite"
    assert not tmp_db.is_file()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus external-alignment ...",
        status="Testing",
        name="Testing missing args",
        method="external-alignment",
        program="",
        version="",
        create_db=True,
    )

    with pytest.raises(
        SystemExit,
        match="ERROR: Missing configuration.extra setting",
    ):
        private_cli.log_external_alignment(database=tmp_db, run_id=1)


def test_bad_config(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check how log-external-alignment handles bad settings."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad.sqlite"
    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus external-alignment ...",
        status="Testing",
        name="Testing bad args",
        method="external-alignment",
        program="",
        version="",
        extra="file=example.fasta;md5=XXX;label=stem",
        create_db=True,
    )

    with pytest.raises(
        SystemExit,
        match="ERROR: configuration.extra='file=example.fasta;md5=XXX;label=stem' unexpected",
    ):
        private_cli.log_external_alignment(database=tmp_db, run_id=1)


def test_missing_alignment(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check how log-external-alignment handles bad settings."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad.sqlite"

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus external-alignment ...",
        status="Testing",
        name="Testing missing alignment file",
        method="external-alignment",
        program="",
        version="",
        extra="md5=XXX;label=stem;alignment=does-not-exist.fasta",
        create_db=True,
    )

    with pytest.raises(
        SystemExit,
        match="ERROR: Missing alignment file .*/does-not-exist.fasta",
    ):
        private_cli.log_external_alignment(database=tmp_db, run_id=1)


def test_bad_checksum(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check how log-external-alignment handles bad settings."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad.sqlite"
    assert not tmp_db.is_file()

    tmp_alignment = tmp_dir / "example.fasta"
    tmp_alignment.touch()

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus external-alignment ...",
        status="Testing",
        name="Testing bad checksum",
        method="external-alignment",
        program="",
        version="",
        extra=f"md5=XXX;label=stem;alignment={tmp_alignment}",
        create_db=True,
    )

    with pytest.raises(
        SystemExit,
        match="ERROR: MD5 checksum of .*/example.fasta didn't match.",
    ):
        private_cli.log_external_alignment(database=tmp_db, run_id=1)


def test_bad_alignment(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Broken alignment as input."""
    tmp_db = Path(tmp_path) / "broken.db"
    assert not tmp_db.is_file()

    tmp_alignment = Path(tmp_path) / "broken.fasta"
    with tmp_alignment.open("w") as handle:
        handle.write("""\
>OP073605 mock 10bp fragment
AACC
>MGV-GENOME-0266457 mock 10bp fragment
AAC
>MGV-GENOME-0264574 mock 9bp fragment
AA
    """)
    md5 = file_md5sum(tmp_alignment)

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus external-alignment ...",
        status="Testing",
        name="Testing bad alignment",
        method="external-alignment",
        program="",
        version="",
        extra=f"md5={md5};label=stem;alignment={tmp_alignment}",
        create_db=True,
    )

    with pytest.raises(
        ValueError,
        match=r"zip\(\) argument 2 is shorter than argument 1",
    ):
        private_cli.log_external_alignment(database=tmp_db, run_id=1)
