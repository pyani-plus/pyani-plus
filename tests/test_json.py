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
"""Tests for the import/export of comparisons in JSON.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pathlib import Path

import pytest

from pyani_plus import db_orm, private_cli, setup_logger


def test_json_import_errors_core(input_genomes_tiny: Path, tmp_path: str) -> None:
    """Confirm expected errors import JSON comparisons."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "json.sqlite"
    assert not tmp_db.is_file()

    tmp_json = tmp_dir / "x.json"
    tmp_json.touch()

    # Missing DB
    with pytest.raises(SystemExit, match=f"Database '{tmp_db}' does not exist"):
        private_cli.import_comparisons(
            tmp_db, json=[tmp_json], debug=False, log=Path("-")
        )

    # Empty DB
    tmp_db.touch()
    with pytest.raises(SystemExit, match="does not contain any configurations"):
        private_cli.import_comparisons(
            tmp_db, json=[tmp_json], debug=False, log=Path("-")
        )

    # Add a config
    private_cli.log_configuration(
        tmp_db,
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        create_db=True,
    )
    with pytest.raises(SystemExit, match="does not contain any genomes"):
        private_cli.import_comparisons(
            tmp_db, json=[tmp_json], debug=False, log=Path("-")
        )

    # Add some genomes
    private_cli.log_genome(
        database=tmp_db,
        fasta=[
            input_genomes_tiny / "MGV-GENOME-0264574.fas",
            input_genomes_tiny / "MGV-GENOME-0266457.fna",
        ],
    )

    # Bad JSON
    with tmp_json.open("w") as handle:
        handle.write("[")
    with pytest.raises(SystemExit, match=f"JSON file '{tmp_json}' invalid"):
        private_cli.import_comparisons(
            tmp_db, json=[tmp_json], debug=False, log=Path("-")
        )

    # Unexpected JSON
    with tmp_json.open("w") as handle:
        handle.write("[]")
    with pytest.raises(
        SystemExit, match=f"JSON file '{tmp_json}' does not use the expected structure"
    ):
        private_cli.import_comparisons(
            tmp_db, json=[tmp_json], debug=False, log=Path("-")
        )


def test_json_import_errors(
    caplog: pytest.LogCaptureFixture, input_genomes_tiny: Path, tmp_path: str
) -> None:
    """Confirm expected errors import JSON comparisons."""
    tmp_dir = Path(tmp_path)
    tmp_db = tmp_dir / "json.sqlite"
    assert not tmp_db.is_file()
    tmp_json = tmp_dir / "x.json"

    # Add a config
    private_cli.log_configuration(
        tmp_db,
        method="guessing",
        program="guestimate",
        version="0.1.2beta3",
        fragsize=100,
        kmersize=51,
        create_db=True,
    )

    # Add some genomes
    private_cli.log_genome(
        database=tmp_db,
        fasta=[
            input_genomes_tiny / "MGV-GENOME-0264574.fas",
            input_genomes_tiny / "MGV-GENOME-0266457.fna",
        ],
    )

    # Empty JSON
    tmp_json.touch()
    caplog.clear()
    private_cli.import_comparisons(tmp_db, json=[tmp_json], debug=True, log=Path("-"))
    output = caplog.text
    assert f"JSON file '{tmp_json}' is empty" in output
    assert f"Imported 0 from '{tmp_json}'" in output

    # Make a mismatched JSON config:
    with tmp_json.open("w") as handle:
        handle.write("""\
{
    "configuration":{"method":"ANIm", "program":"nucmer","version":"3.1", "mode":"mum"},
    "uname":{"system":"Darwin", "release":"24.3.0", "machine":"arm64"},
    "comparisons":[]
}""")

    with pytest.raises(
        SystemExit, match=f"JSON file '{tmp_json}' configuration not in database"
    ):
        private_cli.import_comparisons(
            tmp_db, json=[tmp_json], debug=False, log=Path("-")
        )

    # Partial uname
    with tmp_json.open("w") as handle:
        handle.write("""\
{
    "configuration":{"method":"ANIm", "program":"nucmer","version":"3.1", "mode":"mum"},
    "uname":{"system":"Darwin", "release":"24.3.0"},
    "comparisons":[]
}""")
    with pytest.raises(SystemExit, match=f"JSON file '{tmp_json}' uname incomplete"):
        private_cli.import_comparisons(
            tmp_db, json=[tmp_json], debug=False, log=Path("-")
        )

    # Partial config
    with tmp_json.open("w") as handle:
        handle.write("""\
{
    "configuration":{"method":"ANIm"},
    "uname":{"system":"Darwin", "release":"24.3.0", "machine":"arm64"},
    "comparisons":[]
}""")
    with pytest.raises(
        SystemExit, match=f"JSON file '{tmp_json}' configuration incomplete"
    ):
        private_cli.import_comparisons(
            tmp_db, json=[tmp_json], debug=False, log=Path("-")
        )

    # Match the config, but no comparisons (should be a warning only)
    with tmp_json.open("w") as handle:
        handle.write("""\
{
    "configuration":{"method":"guessing", "program":"guestimate","version":"0.1.2beta3", "fragsize":100,"kmersize":51},
    "uname":{"system":"Darwin", "release":"24.3.0", "machine":"arm64"},
    "comparisons":[]
}""")
    caplog.clear()
    private_cli.import_comparisons(tmp_db, json=[tmp_json], debug=False, log=Path("-"))
    text = caplog.text
    assert f"JSON file '{tmp_json}' has no comparisons" in text

    # Incomplete comparison -- need at least the hashes and identity
    with tmp_json.open("w") as handle:
        handle.write("""\
{
    "configuration":{"method":"guessing", "program":"guestimate","version":"0.1.2beta3", "fragsize":100,"kmersize":51},
    "uname":{"system":"Darwin", "release":"24.3.0", "machine":"arm64"},
    "comparisons":[{"query_hash":"689d3fd6881db36b5e08329cf23cecdd", "identity":0.99}]
}""")
    with pytest.raises(
        SystemExit, match=rf"JSON file '{tmp_json}' comparison\(s\) incomplete"
    ):
        private_cli.import_comparisons(
            tmp_db, json=[tmp_json], debug=False, log=Path("-")
        )

    # Finally, a valid example!
    with tmp_json.open("w") as handle:
        handle.write("""\
{
    "configuration":{"method":"guessing", "program":"guestimate","version":"0.1.2beta3", "fragsize":100,"kmersize":51},
    "uname":{"system":"Darwin", "release":"24.3.0", "machine":"arm64"},
    "comparisons":[
        {
            "query_hash":"689d3fd6881db36b5e08329cf23cecdd",
            "subject_hash":"78975d5144a1cd12e98898d573cf6536",
            "identity":0.99
        }
    ]
}""")
    private_cli.import_comparisons(tmp_db, json=[tmp_json], debug=False, log=Path("-"))

    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, tmp_db)
    assert session.query(db_orm.Comparison).count() == 1
