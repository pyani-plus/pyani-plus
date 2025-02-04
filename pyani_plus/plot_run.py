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
from math import log, nan
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
from rich.progress import Progress

from pyani_plus import PROGRESS_BAR_COLUMNS, db_orm

GREY = (0.7, 0.7, 0.7)
DULL_BLUE = (0.137, 0.412, 0.737)
WHITE = (1.0, 1.0, 1.0)
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

colormaps.register(
    LinearSegmentedColormap.from_list(
        "BuRd",  # blue to red
        (
            (0.0, DULL_BLUE),
            (0.5, WHITE),
            (1.0, DULL_BLUE),
        ),
    )
)


def plot_heatmap(  # noqa: PLR0913
    matrix: pd.DataFrame,
    outdir: Path,
    name: str,  # e.g. tANI
    method: str,
    color_scheme: str,
    formats: tuple[str, ...] = ("tsv", "png", "jpg", "svg", "pdf"),
) -> int:
    """Plot heatmaps for the given matrix."""
    # Can't use square=True with seaborn clustermap, and when clustering
    # asymmetric matrices can end up with different max-length labels
    # for rows vs columns, which distorts layout (non-square heatmap).
    #
    # Decide on figure layout size: a minimum size is required for
    # aesthetics, and a maximum to avoid core dumps on rendering.
    # If we hit the maximum size, we should modify font size.
    maxfigsize = 120
    calcfigsize = matrix.shape[0] * 1.1
    figsize = min(max(8, calcfigsize), maxfigsize)
    if figsize == maxfigsize:
        scale = maxfigsize / calcfigsize
        sns.set_context("notebook", font_scale=scale)

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=(
                "The symmetric non-negative hollow observation matrix"
                " looks suspiciously like an uncondensed distance matrix"
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
            vmax=5 if name == "tANI" else 1,
            figsize=(figsize, figsize),
            linewidths=0.25,
        )

    for ext in formats:
        filename = outdir / f"{method}_{name}_heatmap.{ext}"
        if ext == "tsv":
            # Apply the clustering reordering to match the figure:
            matrix = matrix.iloc[
                figure.dendrogram_row.reordered_ind,
                figure.dendrogram_row.reordered_ind,
            ]
            matrix.to_csv(filename, sep="\t")
        else:
            figure.savefig(filename)
    return len(formats)


def plot_distribution(
    matrix: pd.DataFrame,
    outdir: Path,
    name: str,
    method: str,
    formats: tuple[str, ...] = ("tsv", "png", "jpg", "svg", "pdf"),
) -> int:
    """Plot score distribution and headmap for give matrix.

    Returns the number of plots (should equal number of formats, or zero).
    """
    values = matrix.values.flatten()  # noqa: PD011

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

    # Default rug height=.025 but that can obscure low values.
    # Negative means below the axis instead, needs clip_on=False
    # Adding alpha to try to reveal the density.
    sns.rugplot(
        values,
        ax=axes[1],
        color=rug,
        height=-0.025,
        clip_on=False,
        alpha=0.1,
    )

    # Modify axes after data is plotted
    for _ in axes:
        if name in ["hadamard", "coverage"]:
            _.set_xlim(0, 1.01)
        if name in ["tANI"]:
            _.set_xlim(0, 5.01)
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
    return len(formats)


def plot_scatter(
    run: db_orm.Run,
    outdir: Path,
    formats: tuple[str, ...] = ("tsv", "png", "jpg", "svg", "pdf"),
) -> int:
    """Plot query coverage vs identity for the given run.

    Returns the number of distributions drawn (usually this will be one, zero
    if either property is all nulls).
    """
    method = run.configuration.method
    pairs = [
        (comp.identity, comp.cov_query, comp.aln_length) for comp in run.comparisons()
    ]
    count = len(pairs)
    if not count:
        msg = f"ERROR: No comparisons values in {method} run {run.run_id}"
        sys.exit(msg)
    values = [(x, y, c) for (x, y, c) in pairs if x is not None and y is not None]
    if not values:
        msg = f"WARNING: No valid identity, query-coverage values from {method} run\n"
        sys.stderr.write(msg)
        return 0
    sys.stderr.write(f"Plotting {len(values)}/{count} comparisons\n")

    x_values = [x for (x, y, c) in values]
    y_values = [y for (x, y, c) in values]
    c_values = [c for (x, y, c) in values]

    sys.stderr.write(f"Identity range {min(x_values)} to {max(x_values)}\n")
    sys.stderr.write(f"Query coverage range {min(y_values)} to {max(y_values)}\n")

    # Create the plot
    axes = sns.regplot(
        x=x_values,
        y=y_values,
        fit_reg=False,
        scatter=True,
        scatter_kws={"s": 2, "c": c_values, "color": None},
    )
    axes.set(xlabel="Percent identity (ANI)", ylabel="Query coverage")
    plt.title(f"Percentage identity vs query coverage ({method})")

    for ext in formats:
        filename = outdir / f"{method}_scatter.{ext}"
        if ext == "tsv":
            pass
        else:
            axes.get_figure().savefig(filename)

    return len(formats)


def plot_single_run(  # noqa: C901, PLR0912, PLR0915
    run: db_orm.Run,
    outdir: Path,
    label: str,
    formats: tuple[str, ...] = ("tsv", "png", "jpg", "svg", "pdf"),
) -> int:
    """Plot distributions and heatmaps for given run.

    Draws identity, coverage, hadamard, and tRNA plots of (score distributions
    and heatmaps) for the given run.

    Shows a progress bar in terms of number of scores and plot-types (i.e.
    4 scores times 2 plots giving 8 steps, plus 1 scatter plot).

    Returns number of images drawn.
    """
    method = run.configuration.method
    scores_and_color_schemes = [
        ("identity", "spbnd_BuRd"),
        ("query_cov", "BuRd"),
        ("hadamard", "viridis"),
        ("tANI", "viridis_r"),  # must follow hadamard!
    ]
    did_any_heatmaps = False
    with Progress(*PROGRESS_BAR_COLUMNS) as progress:
        task = progress.add_task(
            "Plotting", total=len(scores_and_color_schemes) * 2 + 1
        )
        # This will print any warnings itself:
        done = plot_scatter(run, outdir)
        progress.advance(task)

        for name, color_scheme in scores_and_color_schemes:
            # The matrices are large, so load them one at a time
            if name == "identity":
                matrix = run.identities
            elif name == "query_cov":
                matrix = run.cov_query
            elif name == "hadamard":
                matrix = run.hadamard

            if matrix is None:
                # This is mainly for mypy to assert the matrix is not None
                msg = f"ERROR: Could not load run {method} {name} matrix"  # pragma: no cover
                sys.exit(msg)  # pragma: no cover

            if name == "tANI":
                # Using run.tani would reload Hadamard and then log transform it.
                # We already loaded the Hadamard matrix, so can transform it here:
                matrix = matrix.map(lambda x: -log(x) if x else nan, na_action="ignore")
            else:
                try:
                    matrix = run.relabelled_matrix(matrix, label)
                except ValueError as err:
                    msg = f"ERROR: {err}"
                    sys.exit(msg)

            nulls = int(matrix.isnull().sum().sum())  # noqa: PD003
            n = len(matrix)
            if nulls == n**2:
                msg = f"WARNING: Cannot plot {name} as all NA\n"
                sys.stderr.write(msg)
                progress.advance(task)  # skipping distribution plots
                progress.advance(task)  # skipping heatmap
                continue

            if matrix.min().min() == matrix.max().max():
                msg = f"WARNING: Skipping {name} plots as all {matrix.min().min()}\n"
                sys.stderr.write(msg)
                progress.advance(task)  # skipping distribution plots
                progress.advance(task)  # skipping heatmap
                continue

            done += plot_distribution(matrix, outdir, name, method, formats)
            progress.advance(task)

            if nulls:
                msg = (
                    f"WARNING: Cannot plot {name} as matrix contains {nulls} nulls"
                    f" (out of {n}²={n**2} {method} comparisons)\n"
                )
                sys.stderr.write(msg)
            else:
                done += plot_heatmap(
                    matrix, outdir, name, method, color_scheme, formats
                )
                did_any_heatmaps = True
            progress.advance(task)
    if not did_any_heatmaps:
        msg = "ERROR: Unable to plot any heatmaps (check for nulls)"
        sys.exit(msg)
    return done
