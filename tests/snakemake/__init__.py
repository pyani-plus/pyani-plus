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
"""Module providing tests of snakemake operation."""

from pathlib import Path

import pandas as pd

from pyani_plus import db_orm, setup_logger


def compare_matrix(
    matrix_df: pd.DataFrame, matrix_path: Path, absolute_tolerance: float | None = None
) -> None:
    """Compare output matrix to expected values from given TSV file.

    The expected output files should now all be using filename stems
    as labels (matching our default export-run output). However, the
    sort order may differ.

    Any absolute_tolerance is only used for floats.
    """

    def strip_colon(text: str) -> str:
        """Drop anything after a colon."""
        assert ":" not in text, text
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
    assert list(matrix_df.columns) == list(expected_df.columns), (
        f"{list(matrix_df.columns)} vs {list(expected_df.columns)}"
    )
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


def label_with_stems(matrix: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    """Relabel and sort a dataframe matrix."""
    return (
        matrix.rename(index=mapping, columns=mapping)
        .sort_index(axis=0)
        .sort_index(axis=1)
    )


def compare_db_matrices(
    database_path: Path,
    matrices_path: Path,
    absolute_tolerance: float = 2e-8,
) -> None:
    """Compare the matrices in the given DB to legacy output from pyANI.

    Assumes there is one and only one run in the database. Checks there
    is one and only one configuration in the database (as a common
    failure was comparisons being logged to a different configuration).

    Assumes the expected matrices on disk are named ``{method}_*.tsv``
    and using MD5 captions internally:

    * ``{method}_aln_lengths.tsv``
    * ``{method}_coverage.tsv``
    * ``{method}_hadamard.tsv``
    * ``{method}_identity.tsv``
    * ``{method}_sim_errors.tsv``

    If any of the files are missing, that comparison is skipped.

    The absolute_tolerance is only used for the floating point matrices.
    """
    logger = setup_logger(None)
    session = db_orm.connect_to_db(logger, database_path)
    run = session.query(db_orm.Run).one()
    assert session.query(db_orm.Configuration).count() == 1, (
        f"Expected one configuration, not {session.query(db_orm.Configuration).count()}"
    )
    n = run.genomes.count()
    assert run.comparisons().count() == n**2, (
        f"Expected {n}²={n**2} comparisons, not {run.comparisons().count()}"
    )
    run.cache_comparisons()  # could be all NaN if constructed in steps
    assert run.identities is not None
    assert matrices_path.is_dir()
    method = run.configuration.method

    md5_to_stem = {_.genome_hash: Path(_.fasta_filename).stem for _ in run.fasta_hashes}

    checked = False
    if (matrices_path / f"{method}_identity.tsv").is_file():
        compare_matrix(
            label_with_stems(run.identities, md5_to_stem),
            matrices_path / f"{method}_identity.tsv",
            absolute_tolerance=absolute_tolerance,
        )
        checked = True
    if (matrices_path / f"{method}_aln_lengths.tsv").is_file():
        compare_matrix(
            label_with_stems(run.aln_length, md5_to_stem),
            matrices_path / f"{method}_aln_lengths.tsv",
        )
        checked = True
    if (matrices_path / f"{method}_coverage.tsv").is_file():
        compare_matrix(
            label_with_stems(run.cov_query, md5_to_stem),
            matrices_path / f"{method}_coverage.tsv",
            absolute_tolerance=absolute_tolerance,
        )
        checked = True
    if (matrices_path / f"{method}_hadamard.tsv").is_file():
        compare_matrix(
            label_with_stems(run.hadamard, md5_to_stem),
            matrices_path / f"{method}_hadamard.tsv",
            absolute_tolerance=absolute_tolerance,
        )
        checked = True
    if (matrices_path / f"{method}_sim_errors.tsv").is_file():
        if method == "dnadiff" and "/viral_example/matrices" in str(matrices_path):
            compare_matrix(
                label_with_stems(run.sim_errors, md5_to_stem),
                matrices_path / f"{method}_sim_errors.tsv",
                absolute_tolerance=1.33,
            )
            checked = True
        else:
            compare_matrix(
                label_with_stems(run.sim_errors, md5_to_stem),
                matrices_path / f"{method}_sim_errors.tsv",
            )
        checked = True
    assert checked, f"Missing expected matrices in {matrices_path}"
