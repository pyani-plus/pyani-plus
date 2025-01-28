# The MIT License
#
# Copyright (c) 2019-2025 University of Strathclyde
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

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap

from pyani_plus import db_orm

GREY = (0.7, 0.7, 0.7)
BLUE = (0.0, 0.0, 1.0)
DULL_BLUE = (0.137, 0.412, 0.737)
WHITE = (1.0, 1.0, 1.0)
RED = (1, 0, 0.0, 0.0)
DULL_RED = (0.659, 0.216, 0.231)

# Custom Matplotlib colourmaps
colormaps.register(
    LinearSegmentedColormap.from_list(
        "spbnd_BuRd",  # species boundary - blue to red
        (
            (0.00, GREY),  # 0% grey
            (0.80, GREY),  # 80% grey (not meaningful)
            (0.80, DULL_BLUE),  # 80% blue
            (0.95, WHITE),  # 95% white (species boundary)
            (1.00, DULL_RED),  # 100% red
        ),
    )
)

# Hadamard is identity * query-coverage, so lower thresholds
colormaps.register(
    LinearSegmentedColormap.from_list(
        "hadamard_BuRd",  # hadamard - blue to red
        (
            (0.00, GREY),  # 0% grey
            (0.25, GREY),  # 25% grey (0.5 * 0.5 = 0.25)
            (0.25, DULL_BLUE),  # 25% blue
            (0.64, WHITE),  # 64% white (0.8 * 0.8 = 0.64)
            (1.00, DULL_RED),  # 100% red
        ),
    )
)

colormaps.register(
    LinearSegmentedColormap.from_list(
        "BuRd",  # blue to red
        (
            (0.0, BLUE),
            (0.5, WHITE),
            (1.0, RED),
        ),
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
            filename = outdir / f"{method}_{name}_heatmap.{ext}"
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


def plot_distributions(
    run: db_orm.Run,
    outdir: Path,
    formats: tuple[str, ...] = ("tsv", "png", "jpg", "svg", "pdf"),
) -> int:
    """Plot identity, coverage, and hadamard distributions for the given run.

    Returns the number of distributions drawn (skips any which are all nulls).
    """
    method = run.configuration.method
    plots_done = 0
    for matrix, name, color_scheme in (
        (run.identities, "identity", "spbnd_BuRd"),
        (run.cov_query, "query_cov", "BuRd"),
        (run.hadamard, "hadamard", "hadamard_BuRd"),
    ):
        del color_scheme  # should use this on the left plot?

        if matrix is None:
            # This is mainly for mypy to assert the matrix is not None
            msg = f"ERROR: Could not load run {method} matrix"  # pragma: no cover
            sys.exit(msg)  # pragma: no cover

        nulls = int(matrix.isnull().sum().sum())
        n = len(matrix)
        if nulls == n**2:
            msg = f"WARNING: Cannot plot {name} as all NA\n"
            sys.stderr.write(msg)
            continue

        values = matrix.values.flatten()

        fill = "#A6C8E0"
        rug = "#2678B2"
        figure, axes = plt.subplots(1, 2, figsize=(15, 5))
        figure.suptitle(f"{name} distribution")
        sns.histplot(
            values,
            ax=axes[0],
            stat="count",
            element="step",
            color=fill,
            edgecolor=fill,
        )
        axes[0].set_ylim(ymin=0)
        sns.kdeplot(values, ax=axes[1], warn_singular=False)
        sns.rugplot(values, ax=axes[1], color=rug)

        # Modify axes after data is plotted
        for _ in axes:
            if name in ["hadamard", "coverage"]:
                _.set_xlim(0, 1.01)
            elif name == "identity":
                # 80% default matches the heatmap grey/blue default
                _.set_xlim(0.80, 1.01)

        # Tidy figure
        figure.tight_layout(rect=[0, 0.03, 1, 0.95])

        for ext in formats:
            filename = outdir / f"{method}_{name}_dist.{ext}"
            if ext == "tsv":
                pass
            else:
                figure.savefig(filename)
        plots_done += 1
    return plots_done
