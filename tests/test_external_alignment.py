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

import logging
from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, public_cli, setup_logger
from pyani_plus.utils import file_md5sum

# Listing the fragments in MD5 order to match the matrix in DB:
MOCK_3_BY_11_ALIGNMENT = """\
>OP073605 mock 10bp fragment for 5584c7029328dc48d33f95f0a78f7e57
GACC-GGTTTT
>MGV-GENOME-0264574 mock 9bp fragment for 689d3fd6881db36b5e08329cf23cecdd
AACC-GG-TTT
>MGV-GENOME-0266457 mock 10bp fragment for 78975d5144a1cd12e98898d573cf6536
AACC-GGATTT
"""
# Consider this pair (1st and 2nd entries):
#
# 5584c7029328dc48d33f95f0a78f7e57 GACC-GGTTTT query length 10
# 689d3fd6881db36b5e08329cf23cecdd AACC-GG-TTT subject length 9
#
# Equivalent to this:
#
# 5584c7029328dc48d33f95f0a78f7e57 GACCGGTTTT query length 10
# 689d3fd6881db36b5e08329cf23cecdd AACCGG-TTT subject length 9
#
# 8 matches 1 mismatch in alignment of gapped-length 10, so identity 8/10=0.8
# However query coverage (8+1)/10 = 0.9, subject cover (8+1)/9 = 1.0
#
# So identity [1.0 0.8 ...]    query coverage [1.0 0.9 ...]
#             [0.8 1.0 ...]                   [1.0 1.0 ...]
#             [... ... 1.0]                   [... ... 1.0]
#
# Now consider this pair (1st and 3rd entries):
#
# 5584c7029328dc48d33f95f0a78f7e57 GACC-GGTTTT query length 10
# 78975d5144a1cd12e98898d573cf6536 AACC-GGATTT subject length 10
#
# Equivalent to this:
#
# 5584c7029328dc48d33f95f0a78f7e57 GACCGGTTTT
# 78975d5144a1cd12e98898d573cf6536 AACCGGATTT
#
# 8 matches 2 mismatch in an alignment of length 10, identity 8/10=0.8
# Both query and subject full coverage (8+2)/10, so coverage 1.0
#
# So identity [1.0 ... 0.8]    query coverage [1.0 ... 1.0]
#             [... 1.0 ...]                   [... 1.0 ...]
#             [0.8 ... 1.0]                   [1.0 ... 1.0]
#
# Finally, consider this pair (2nd and 3rd entries):
#
# 689d3fd6881db36b5e08329cf23cecdd AACC-GG-TTT query length 9
# 78975d5144a1cd12e98898d573cf6536 AACC-GGATTT subject length 10
#
# Equivalent to this:
#
# 689d3fd6881db36b5e08329cf23cecdd AACCGG-TTT query length 9
# 78975d5144a1cd12e98898d573cf6536 AACCGGATTT subject length 10
#
# 9 matches 0 mismatch in alignment of length 10, so identity 9/10 = 0.9
# However query coverage (9+0)/9 = 1.0, subject cover (9+0)/10 = 0.9
#
# So identity [1.0 ... ...]    query coverage [1.0 ... ...]
#             [... 1.0 0.9]                   [... 1.0 1.0]
#             [... 0.9 1.0]                   [... 0.9 1.0]
#
# Overall:
#
# So identity [1.0 0.8 0.8]    query coverage [1.0 0.9 1.0]
#             [0.8 1.0 0.9]                   [1.0 1.0 1.0]
#             [0.8 0.9 1.0]                   [1.0 0.9 1.0]
#
MOCK_3_BY_11_DF_IDENTITY = (
    '{"columns":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
    '"index":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
    '"data":[[1.0,0.8,0.8],[0.8,1.0,0.9],[0.8,0.9,1.0]]}'
)
MOCK_3_BY_11_DF_COV_QUERY = (
    '{"columns":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
    '"index":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
    '"data":[[1.0,0.9,1.0],[1.0,1.0,1.0],[1.0,0.9,1.0]]}'
)


def test_simple_mock_alignment_stem(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Mock alignments using filename stem naming.

    These are 9 or 10bp, giving 10bp long pairwise alignments.

    There is a deliberate all-gap column in the alignment to ensure any matching
    gaps are dropped in the pairwise comparisons. So none are 11bp alignments
    despite the alignment having 11 columns.
    """
    tmp_db = Path(tmp_path) / "stems.db"
    assert not tmp_db.is_file()

    tmp_alignment = Path(tmp_path) / "stems.fasta"

    with tmp_alignment.open("w") as handle:
        handle.write(MOCK_3_BY_11_ALIGNMENT)

    caplog.set_level(logging.INFO)
    public_cli.external_alignment(
        input_genomes_tiny, tmp_db, create_db=True, alignment=tmp_alignment
    )
    output = caplog.text
    assert "external-alignment run setup with 3 genomes" in output, output

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    assert session.query(db_orm.Comparison).count() == 9  # noqa: PLR2004

    run = db_orm.load_run(session, 1, check_complete=True)
    assert run.df_identity == MOCK_3_BY_11_DF_IDENTITY
    assert run.df_cov_query == MOCK_3_BY_11_DF_COV_QUERY


def test_simple_mock_alignment_md5(
    caplog: pytest.LogCaptureFixture,
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

    caplog.set_level(logging.INFO)
    public_cli.external_alignment(
        input_genomes_tiny, tmp_db, create_db=True, alignment=tmp_alignment, label="md5"
    )
    output = caplog.text
    assert "external-alignment run setup with 3 genomes" in output, output

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    assert session.query(db_orm.Comparison).count() == 9  # noqa: PLR2004
    run = db_orm.load_run(session, 1, check_complete=True)

    assert run.df_identity == (
        '{"columns":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"index":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"data":[[1.0,0.8,0.8],[0.8,1.0,0.8],[0.8,0.8,1.0]]}'
    )


def test_simple_mock_alignment_filename(
    caplog: pytest.LogCaptureFixture,
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

    caplog.set_level(logging.INFO)
    public_cli.external_alignment(
        input_genomes_tiny,
        tmp_db,
        create_db=True,
        alignment=tmp_alignment,
        label="filename",
    )
    output = caplog.text
    assert "external-alignment run setup with 3 genomes" in output, output

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    assert session.query(db_orm.Comparison).count() == 9  # noqa: PLR2004
    run = db_orm.load_run(session, 1, check_complete=True)
    assert run.df_identity == (
        '{"columns":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"index":["5584c7029328dc48d33f95f0a78f7e57","689d3fd6881db36b5e08329cf23cecdd","78975d5144a1cd12e98898d573cf6536"],'
        '"data":[[1.0,0.8,0.8],[0.8,1.0,0.8],[0.8,0.8,1.0]]}'
    )


def test_resume(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Resume importing an external alignment (from nothing)."""
    tmp_db = Path(tmp_path) / "stems.db"
    assert not tmp_db.is_file()

    tmp_alignment = Path(tmp_path) / "stems.fasta"
    with tmp_alignment.open("w") as handle:
        handle.write(MOCK_3_BY_11_ALIGNMENT)

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

    caplog.set_level(logging.INFO)
    public_cli.resume(tmp_db)
    output = caplog.text
    assert (
        "Database already has 0 of 3²=9 external-alignment comparisons, 9 needed"
        in output
    ), output

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    assert session.query(db_orm.Comparison).count() == 9  # noqa: PLR2004
    run = db_orm.load_run(session, 1, check_complete=True)
    assert run.df_identity == MOCK_3_BY_11_DF_IDENTITY
    assert run.df_cov_query == MOCK_3_BY_11_DF_COV_QUERY


def test_resume_partial(
    caplog: pytest.LogCaptureFixture,
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Resume importing an external alignment (from nothing)."""
    tmp_db = Path(tmp_path) / "stems.db"
    assert not tmp_db.is_file()

    tmp_alignment = Path(tmp_path) / "stems.fasta"
    with tmp_alignment.open("w") as handle:
        handle.write(MOCK_3_BY_11_ALIGNMENT)

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
    private_cli.log_comparison(
        database=tmp_db,
        config_id=1,
        query_fasta=input_genomes_tiny / "OP073605.fasta",
        subject_fasta=input_genomes_tiny / "OP073605.fasta",
        identity=1.0,
        aln_length=10,
        sim_errors=0,
        cov_query=1.0,
        cov_subject=1.0,
    )  # values as per test_simple_mock_alignment_stem

    caplog.set_level(logging.INFO)
    public_cli.resume(tmp_db)
    output = caplog.text
    assert (
        "Database already has 1 of 3²=9 external-alignment comparisons, 8 needed"
        in output
    ), output

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    assert session.query(db_orm.Comparison).count() == 9  # noqa: PLR2004
    run = db_orm.load_run(session, 1, check_complete=True)
    assert run.df_identity == MOCK_3_BY_11_DF_IDENTITY
    assert run.df_cov_query == MOCK_3_BY_11_DF_COV_QUERY


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
        match="We expect no tool information, but run-id 1 used should-be-blank version 1.0 instead.",
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
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "wrong-method.sqlite"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "x.json"

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

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.load_run(session, 1)

    with pytest.raises(
        SystemExit,
        match="Run-id 1 expected guessing results",
    ):
        private_cli.compute_external_alignment(
            logger, tmp_dir, session, run, tmp_json, tmp_dir, {}, {}, {}, ""
        )
    session.close()


def test_bad_program(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check log-external-alignment rejects a run for another method."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "wrong-method.sqlite"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "bad.json"

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus external-alignment ...",
        status="Testing",
        name="Testing external-alignment",
        method="external-alignment",
        program="should-be-blank",
        version="",
        create_db=True,
    )

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.load_run(session, 1)

    with pytest.raises(
        SystemExit,
        match="configuration.program='should-be-blank' unexpected",
    ):
        private_cli.compute_external_alignment(
            logger, tmp_dir, session, run, tmp_json, tmp_dir, {}, {}, {}, ""
        )
    session.close()


def test_bad_version(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check log-external-alignment rejects a run for another method."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "wrong-method.sqlite"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "bad.json"

    private_cli.log_run(
        fasta=input_genomes_tiny,
        database=tmp_db,
        cmdline="pyani-plus external-alignment ...",
        status="Testing",
        name="Testing external-alignment",
        method="external-alignment",
        program="",
        version="should-be-blank",
        create_db=True,
    )

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.load_run(session, 1)

    with pytest.raises(
        SystemExit,
        match="configuration.version='should-be-blank' unexpected",
    ):
        private_cli.compute_external_alignment(
            logger, tmp_dir, session, run, tmp_json, tmp_dir, {}, {}, {}, ""
        )
    session.close()


def test_no_config(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check how log-external-alignment handles bad settings."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad.sqlite"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "no.json"

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

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.load_run(session, 1)

    with pytest.raises(
        SystemExit,
        match="Missing configuration.extra setting",
    ):
        private_cli.compute_external_alignment(
            logger, tmp_dir, session, run, tmp_json, tmp_dir, {}, {}, {}, ""
        )
    session.close()


def test_bad_config(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check how log-external-alignment handles bad settings."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad.sqlite"
    tmp_json = tmp_dir / "bad.json"
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

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.load_run(session)

    with pytest.raises(
        SystemExit,
        match="configuration.extra='file=example.fasta;md5=XXX;label=stem' unexpected",
    ):
        private_cli.compute_external_alignment(
            logger, tmp_dir, session, run, tmp_json, tmp_dir, {}, {}, {}, ""
        )
    session.close()


def test_missing_alignment(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check how log-external-alignment handles bad settings."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad.sqlite"
    tmp_json = tmp_dir / "bad.json"

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

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.load_run(session)

    with pytest.raises(
        SystemExit,
        match="Missing alignment file .*/does-not-exist.fasta",
    ):
        private_cli.compute_external_alignment(
            logger, tmp_dir, session, run, tmp_json, tmp_dir, {}, {}, {}, ""
        )
    session.close()


def test_bad_checksum(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Check how log-external-alignment handles bad settings."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "bad.sqlite"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "bad.json"

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

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.load_run(session)

    with pytest.raises(
        SystemExit,
        match="MD5 checksum of .*/example.fasta didn't match.",
    ):
        private_cli.compute_external_alignment(
            logger, tmp_dir, session, run, tmp_json, tmp_dir, {}, {}, {}, ""
        )
    session.close()


def test_bad_alignment(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """Broken alignment as input."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "broken.db"
    assert not tmp_db.is_file()
    json_file = tmp_dir / "bad.json"

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

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.load_run(session)

    with pytest.raises(
        SystemExit,
        match=(
            "Bad external-alignment, different lengths 3 and 2"
            " from MGV-GENOME-0266457 and MGV-GENOME-0264574"
        ),
    ):
        private_cli.compute_external_alignment(
            logger,
            tmp_dir,
            session,
            run,
            json_file,
            tmp_dir,
            {},
            {},
            {
                "5584c7029328dc48d33f95f0a78f7e57": 57793,
                "689d3fd6881db36b5e08329cf23cecdd": 39253,
                "78975d5144a1cd12e98898d573cf6536": 39594,
            },
            "689d3fd6881db36b5e08329cf23cecdd",
        )
    session.close()


def test_partial_alignment(
    tmp_path: str,
    input_genomes_tiny: Path,
) -> None:
    """MSA input without all expected genomes."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "broken.db"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "msa.json"

    tmp_alignment = Path(tmp_path) / "broken.fasta"
    with tmp_alignment.open("w") as handle:
        handle.write("""\
>OP073605 mock 10bp fragment
AACC
>MGV-GENOME-0266457 mock 10bp fragment
AACT
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

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    run = db_orm.load_run(session)

    # This one should work...
    private_cli.compute_external_alignment(
        logger,
        tmp_dir,
        session,
        run,
        tmp_json,
        tmp_dir,
        {},
        {},
        {
            "5584c7029328dc48d33f95f0a78f7e57": 57793,
            "689d3fd6881db36b5e08329cf23cecdd": 39253,
            "78975d5144a1cd12e98898d573cf6536": 39594,
        },
        "5584c7029328dc48d33f95f0a78f7e57",
    )
    with pytest.raises(
        SystemExit,
        match="Did not find subject 689d3fd6881db36b5e08329cf23cecdd in broken.fasta",
    ):
        private_cli.compute_external_alignment(
            logger,
            tmp_dir,
            session,
            run,
            tmp_json,
            tmp_dir,
            {},
            {},
            {
                "5584c7029328dc48d33f95f0a78f7e57": 57793,
                "689d3fd6881db36b5e08329cf23cecdd": 39253,
                "78975d5144a1cd12e98898d573cf6536": 39594,
            },
            "689d3fd6881db36b5e08329cf23cecdd",
        )
    session.close()
