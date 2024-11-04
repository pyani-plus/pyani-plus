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
"""Module providing tests of snakemake operation."""

from pathlib import Path

import pandas as pd

from pyani_plus import db_orm


def compare_matrix(
    matrix_df: pd.DataFrame, matrix_path: Path, absolute_tolerance: float | None = None
) -> None:
    """Compare output matrix to expected values from given TSV file.

    The output from legacy pyANI v0.3 should be using MD5 captions,
    but will have appended colon-index to them. Also the order will
    be different, the current ORM returns the matrix sorted by MD5.

    Any absolute_tolerance is only used for floats.
    """

    def strip_colon(text: str) -> str:
        """Drop anything after a colon."""
        return text.split(":", 1)[0] if ":" in text else text

    # Using converters to fix the row names (in column 0)
    # and rename method to fix the column names:
    expected_df = (
        pd.read_csv(
            matrix_path, sep="\t", header=0, index_col=0, converters={0: strip_colon}
        )
        .rename(columns=strip_colon)
        .sort_index(axis=0)
        .sort_index(axis=1)
    )
    assert list(matrix_df.columns) == list(expected_df.columns)
    if not expected_df.dtypes.equals(matrix_df.dtypes):
        # This happens with some old pyANI output using floats for ints
        # Cast both to float
        expected_df = expected_df.astype(float)
        matrix_df = matrix_df.astype(float)
    if absolute_tolerance is None:
        pd.testing.assert_frame_equal(matrix_df, expected_df, obj=matrix_path.stem)
    else:
        pd.testing.assert_frame_equal(
            matrix_df, expected_df, obj=matrix_path.stem, atol=absolute_tolerance
        )


def compare_matrices(
    database_path: Path, matrices_path: Path, absolute_tolerance: float = 2e-8
) -> None:
    """Compare the matrices in the given DB to legacy output from pyANI.

    Assumes there is one and only one run in the database. Checks there
    is one and only one configuration in the database (as a common
    failure was comparisons being logged to a different configuration).

    Assumes the legacy matrices are named ``matrix_*.tsv`` and use MD5
    captions internally:

    * ``matrix_aln_lengths.tsv``
    * ``matrix_coverage.tsv``
    * ``matrix_hadamard.tsv``
    * ``matrix_identity.tsv``
    * ``matrix_sim_errors.tsv``

    If any of the files are missing, the comparison is skipped.

    The absolute_tolerance is only used for the floating point matrices.
    """
    session = db_orm.connect_to_db(database_path)
    run = session.query(db_orm.Run).one()
    assert (
        session.query(db_orm.Configuration).count() == 1
    ), f"Expected one configuration, not {session.query(db_orm.Configuration).count()}"
    n = run.genomes.count()
    assert (
        run.comparisons().count() == n**2
    ), f"Expected {n}x{n}={n**2} comparisons, not {run.comparisons().count()}"

    assert run.identities is not None
    assert matrices_path.is_dir()

    checked = False
    if (matrices_path / "matrix_identity.tsv").is_file():
        compare_matrix(
            run.identities,
            matrices_path / "matrix_identity.tsv",
            absolute_tolerance=absolute_tolerance,
        )
        checked = True
    if (matrices_path / "matrix_aln_lengths.tsv").is_file():
        compare_matrix(run.aln_length, matrices_path / "matrix_aln_lengths.tsv")
        checked = True
    if (matrices_path / "matrix_coverage.tsv").is_file():
        compare_matrix(
            run.cov_query,
            matrices_path / "matrix_coverage.tsv",
            absolute_tolerance=absolute_tolerance,
        )
        checked = True
    if (matrices_path / "matrix_hadamard.tsv").is_file():
        compare_matrix(
            run.hadamard,
            matrices_path / "matrix_hadamard.tsv",
            absolute_tolerance=absolute_tolerance,
        )
        checked = True
    if (matrices_path / "matrix_sim_errors.tsv").is_file():
        if "dnadiff" in str(matrices_path):
            compare_matrix(
                run.sim_errors,
                matrices_path / "matrix_sim_errors.tsv",
                absolute_tolerance=1.33,
            )
            checked = True
        else:
            compare_matrix(run.sim_errors, matrices_path / "matrix_sim_errors.tsv")
        checked = True
    assert checked, f"Missing expected matrices in {matrices_path}"
