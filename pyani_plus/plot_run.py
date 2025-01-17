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
"""Code for plotting a single run (heatmaps etc)."""

import sys
import warnings
from pathlib import Path

import seaborn as sns
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap

from pyani_plus import db_orm

# Custom Matplotlib colourmaps
# 1a) Map for species boundaries (95%: 0.95), blue for values at
# 0.9 or below, red for values at 1.0; white at 0.95.
# Also, anything below 0.7 is 70% grey
colormaps.register(
    LinearSegmentedColormap(
        "spbnd_BuRd",
        {
            "red": (
                (0.0, 0.0, 0.7),
                (0.7, 0.7, 0.0),
                (0.9, 0.0, 0.0),
                (0.95, 1.0, 1.0),
                (1.0, 1.0, 1.0),
            ),
            "green": (
                (0.0, 0.0, 0.7),
                (0.7, 0.7, 0.0),
                (0.9, 0.0, 0.0),
                (0.95, 1.0, 1.0),
                (1.0, 0.0, 0.0),
            ),
            "blue": (
                (0.0, 0.0, 0.7),
                (0.7, 0.7, 1.0),
                # skips (0.9, 0.0, 0.0), in legacy code (why?)
                (0.95, 1.0, 1.0),
                (1.0, 0.0, 0.0),
            ),
        },
    )
)


# 1b) Map for species boundaries (95%: 0.95), blue for values at
# 0.64 (0.8 * 0.8) or below, red for values at 1.0; white at 0.9.
# Also, anything below 0.25 (0.5 * 0.5) is 70% grey
colormaps.register(
    LinearSegmentedColormap(
        "hadamard_BuRd",
        {
            "red": (
                (0.0, 0.0, 0.7),
                (0.25, 0.7, 0.0),
                (0.64, 0.0, 0.0),
                (0.64, 1.0, 1.0),
                (1.0, 1.0, 1.0),
            ),
            "green": (
                (0.0, 0.0, 0.7),
                (0.25, 0.7, 0.0),
                (0.64, 0.0, 0.0),
                (0.64, 1.0, 1.0),
                (1.0, 0.0, 0.0),
            ),
            "blue": (
                (0.0, 0.0, 0.7),
                (0.25, 0.7, 1.0),
                # skips (0.64, 0.0, 0.0), in legacy code (why?)
                (0.64, 1.0, 1.0),
                (1.0, 0.0, 0.0),
            ),
        },
    )
)


# 2) Blue for values at 0.0, red for values at 1.0; white at 0.5
colormaps.register(
    LinearSegmentedColormap(
        "BuRd",
        {
            "red": ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)),
            "green": ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)),
            "blue": ((0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)),
        },
    )
)


def plot_heatmaps(
    run: db_orm.Run,
    outdir: Path,
    label: str,
    formats: tuple[str, ...] = ("tsv", "png", "jpg", "svg", "pdf"),
) -> int:
    """Plot the heatmaps for the given run.

    Returns the number of heatmaps drawn (depends on nulls).
    """
    # Can't use square=True with seaborn clustermap, and when clustering
    # asymmetric matrices can end up with different max-length labels
    # for rows vs columns, which distorts layout (non-square heatmap).
    #
    # Decide on figure layout size: a minimum size is required for
    # aesthetics, and a maximum to avoid core dumps on rendering.
    # If we hit the maximum size, we should modify font size.
    method = run.configuration.method

    maxfigsize = 120
    if run.identities is None:
        # This is mainly for mypy to assert the matrix is not None
        msg = "ERROR: Could not load run identities matrix"  # pragma: no cover
        sys.exit(msg)  # pragma: no cover
    calcfigsize = run.identities.shape[0] * 1.1
    figsize = min(max(8, calcfigsize), maxfigsize)
    if figsize == maxfigsize:
        scale = maxfigsize / calcfigsize
        sns.set_context("notebook", font_scale=scale)

    heatmaps_done = 0
    for matrix, name, color_scheme in (
        (run.identities, "identity", "spbnd_BuRd"),
        (run.cov_query, "query_cov", "BuRd"),
        (run.hadamard, "hadamard", "hadamard_BuRd"),
    ):
        if matrix is None:
            # This is mainly for mypy to assert the matrix is not None
            msg = f"ERROR: Could not load run {method} matrix"  # pragma: no cover
            sys.exit(msg)  # pragma: no cover

        nulls = int(matrix.isnull().sum().sum())
        n = len(matrix)
        if nulls:
            msg = (
                f"WARNING: Cannot plot {name} as matrix contains {nulls} nulls"
                f" (out of {n}²={n**2} {method} comparisons)\n"
            )
            sys.stderr.write(msg)
            continue

        try:
            matrix = run.relabelled_matrix(matrix, label)  # noqa: PLW2901
        except ValueError as err:
            msg = f"ERROR: {err}"
            sys.exit(msg)

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=(
                    "scipy.cluster: The symmetric non-negative hollow observation"
                    " matrix looks suspiciously like an uncondensed distance matrix"
                ),
            )
            warnings.filterwarnings(
                "ignore",
                message=(
                    "Clustering large matrix with scipy. Installing"
                    " `fastcluster` may give better performance."
                ),
            )
            figure = sns.clustermap(
                matrix,
                cmap=color_scheme,
                vmin=0,
                vmax=1,
                figsize=(figsize, figsize),
                linewidths=0.25,
            )

        for ext in formats:
            filename = outdir / f"{method}_{name}.{ext}"
            if ext == "tsv":
                # Apply the clustering reordering to match the figure:
                matrix = matrix.iloc[  # noqa: PLW2901
                    figure.dendrogram_row.reordered_ind,
                    figure.dendrogram_row.reordered_ind,
                ]
                matrix.to_csv(filename, sep="\t")
            else:
                figure.savefig(filename)
        heatmaps_done += 1
    return heatmaps_done
