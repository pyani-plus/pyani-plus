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

import logging
import warnings
from math import ceil, log, nan, sqrt
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib import cm, colormaps, colors
from matplotlib.colors import LinearSegmentedColormap
from rich.progress import Progress
from sqlalchemy.orm import Session

from pyani_plus import GRAPHICS_FORMATS, PROGRESS_BAR_COLUMNS, db_orm, log_sys_exit

mpl.use("agg")  # non-interactive backend

ORANGE = (0.934, 0.422, 0)
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
            (1.0, DULL_RED),
        ),
    )
)


def plot_heatmap(  # noqa: PLR0913
    matrix: pd.DataFrame,
    outdir: Path,
    name: str,  # e.g. tANI
    method: str,
    color_scheme: str,
    formats: tuple[str, ...] = GRAPHICS_FORMATS,
    na_fill: float = 0,
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
    if figsize == maxfigsize:  # pragma: nocover
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
            matrix.fillna(na_fill),
            mask=matrix.isna(),
            cmap=colormaps[color_scheme].with_extremes(bad=ORANGE),
            vmin=0,
            vmax=5 if name == "tANI" else 1,
            figsize=(figsize, figsize),
            linewidths=0.25,
        )

    import sys

    sys.stderr.write(f"{figure.ax_col_dendrogram}\n")
    sys.stderr.write(f"{figure.ax_row_dendrogram}\n")

    # adjust cbar to avoid overlapping with the dendrogram
    row_dendrogram_box = figure.ax_row_dendrogram.get_position()
    col_dendrogram_box = figure.ax_col_dendrogram.get_position()
    figure.ax_cbar.set_position(
        (
            row_dendrogram_box.xmin,
            col_dendrogram_box.ymin,
            min(0.05, row_dendrogram_box.width),
            col_dendrogram_box.height,
        )
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
    values: list[float],
    outdir: Path,
    name: str,
    method: str,
    formats: tuple[str, ...] = GRAPHICS_FORMATS,
) -> int:
    """Plot score distribution for give matrix.

    Returns the number of plots (should equal number of formats, or zero).
    """
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

    # In a possible bug, the rug-plot ignores the x-limits,
    # so mask the data here
    for _ in axes:
        if name in ["hadamard", "coverage"]:
            _.set_xlim(0, 1.01)
            values = [v for v in values if v is not None and 0 <= v <= 1.01]  # noqa: PLR2004
        if name in ["tANI"]:
            _.set_xlim(0, 5.01)
            values = [v for v in values if v is not None and 0 <= v <= 5.01]  # noqa: PLR2004
        elif name == "identity":
            # 80% default matches the heatmap grey/blue default
            _.set_xlim(0.80, 1.01)
            values = [v for v in values if v is not None and 0.80 <= v <= 1.01]  # noqa: PLR2004

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

    # Tidy figure
    figure.tight_layout(rect=[0, 0.03, 1, 0.95])

    for ext in formats:
        filename = outdir / f"{method}_{name}_dist.{ext}"
        if ext == "tsv":
            pass
        else:
            figure.savefig(filename)
    plt.close()
    return len(formats)


def plot_scatter(
    logger: logging.Logger,
    run: db_orm.Run,
    outdir: Path,
    formats: tuple[str, ...] = GRAPHICS_FORMATS,
) -> int:
    """Plot query coverage & tANI vs identity for the given run.

    Returns the number of distributions drawn (usually this will be one per format,
    but zero if either property is all nulls).
    """
    method = run.configuration.method
    # Do query coverage first - if that fails, would not be able to calculate tANI
    for y_caption in ("Query coverage", "tANI"):
        if y_caption == "tANI":
            pairs = [
                (
                    comp.identity,
                    None
                    if comp.identity is None or comp.cov_query is None
                    else -log(comp.identity * comp.cov_query),
                    comp.query.length,
                )
                for comp in run.comparisons()
            ]
        else:
            pairs = [
                (comp.identity, comp.cov_query, comp.query.length)
                for comp in run.comparisons()
            ]
        count = len(pairs)  # including missing values
        values = [(x, y, c) for (x, y, c) in pairs if x is not None and y is not None]
        if not values:
            msg = f"No valid identity, {y_caption} values from {method} run"
            logger.warning(msg)
            return 0
        msg = f"Plotting {len(values)}/{count} {y_caption} vs identity {method} comparisons"
        logger.info(msg)

        x_values = [x for (x, y, c) in values]
        y_values = [y for (x, y, c) in values]
        c_values = [c for (x, y, c) in values]

        msg = f"Identity range {min(x_values)} to {max(x_values)}"
        logger.debug(msg)
        msg = f"{y_caption} range {min(y_values)} to {max(y_values)}"
        logger.debug(msg)

        # Create the plot
        joint_grid = sns.jointplot(
            x=x_values,
            y=y_values,
            kind="scatter",
            joint_kws={"s": 2, "c": c_values, "color": None},
        )
        joint_grid.set_axis_labels(xlabel="Percent identity (ANI)", ylabel=y_caption)
        # Can use plt.title(...) but would need to arrange all the placements

        # Shrink for color bar (aka cbar) on right
        plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)
        plt.colorbar(
            cm.ScalarMappable(norm=colors.Normalize(min(c_values), max(c_values))),
            cax=joint_grid.fig.add_axes([0.85, 0.25, 0.05, 0.4]),
            label="Query length (bp)",
        )

        if y_caption == "Query coverage":
            # avoid spaces in filename:
            y_caption = "query_cov"  # noqa: PLW2901
        for ext in formats:
            filename = outdir / f"{method}_{y_caption}_scatter.{ext}"
            if ext == "tsv":
                with filename.open("w") as handle:
                    handle.write(f"#identity\t{y_caption}\tquery_length\n")
                    for x, y, c in values:
                        handle.write(f"{x}\t{y}\t{c}\n")
            else:
                joint_grid.savefig(filename)
        # Clear plot to avoid over-plotting next image
        plt.close()

    return len(formats)


def plot_single_run(
    logger: logging.Logger,
    run: db_orm.Run,
    outdir: Path,
    label: str,
    formats: tuple[str, ...] = GRAPHICS_FORMATS,
) -> int:
    """Plot distributions and heatmaps for given run.

    Draws identity, coverage, hadamard, and tRNA plots of (score distributions
    and heatmaps) for the given run.

    Shows a progress bar in terms of number of scores and plot-types (i.e.
    4 scores times 2 plots giving 8 steps, plus 2 scatter plots).

    Returns number of images drawn.
    """
    method = run.configuration.method
    scores_and_color_schemes = [
        ("identity", "spbnd_BuRd", 0),
        ("query_cov", "BuRd", 0),
        ("hadamard", "viridis", 0),
        ("tANI", "viridis_r", -5),  # must follow hadamard!
    ]
    with Progress(*PROGRESS_BAR_COLUMNS) as progress:
        task = progress.add_task(
            "Plotting", total=len(scores_and_color_schemes) * 2 + 2
        )
        # This will print any warnings itself:
        done = plot_scatter(logger, run, outdir)
        # Need finer grained progress logging:
        progress.advance(task)
        progress.advance(task)

        for name, color_scheme, na_fill in scores_and_color_schemes:
            # The matrices are large, so load them one at a time
            if name == "identity":
                matrix = run.identities
            elif name == "query_cov":
                matrix = run.cov_query
            elif name == "hadamard":
                matrix = run.hadamard

            if matrix is None:  # pragma: no cover
                # This is mainly for mypy to assert the matrix is not None
                msg = f"Could not load run {method} {name} matrix"
                log_sys_exit(logger, msg)
                return 0  # won't be called but mypy doesn't understand (yet)

            if name == "tANI":
                # Using run.tani would reload Hadamard and then log transform it.
                # We already loaded the Hadamard matrix, so can transform it here:
                matrix = matrix.map(lambda x: -log(x) if x else nan, na_action="ignore")
            else:
                try:
                    matrix = run.relabelled_matrix(matrix, label)
                except ValueError as err:
                    msg = str(err)
                    log_sys_exit(logger, msg)

            nulls = int(matrix.isnull().sum().sum())  # noqa: PD003
            n = len(matrix)
            if nulls == n**2:
                msg = f"Cannot plot {name} as all NA"
                logger.warning(msg)
                progress.advance(task)  # skipping distribution plots
                progress.advance(task)  # skipping heatmap
                continue
            if nulls:
                msg = (
                    f"{name} matrix contains {nulls} nulls"
                    f" (out of {n}²={n**2} {method} comparisons)"
                )
                logger.warning(msg)

            values = matrix.values.flatten()  # noqa: PD011
            done += plot_distribution(values, outdir, name, method, formats)
            del values
            progress.advance(task)

            done += plot_heatmap(
                matrix, outdir, name, method, color_scheme, formats, na_fill
            )
            progress.advance(task)
    return done


def plot_run_comparison(  # noqa: C901, PLR0912, PLR0913, PLR0915
    logger: logging.Logger,
    session: Session,
    run: db_orm.Run,
    other_runs: list[int],
    outdir: Path,
    columns: int = 0,
    field: str = "identity",
    formats: tuple[str, ...] = GRAPHICS_FORMATS,
    hist_bins: int = 30,
) -> int:
    """Plot some identity comparisons between runs."""
    # Offer to do this for other properties, e.g. tANI?

    # Ignore E711, must be != None in SQLalchemy
    reference_values_by_hash = {
        (_.query_hash, _.subject_hash): _.identity
        for _ in run.comparisons().where(db_orm.Comparison.identity != None)  # noqa: E711
    }
    queries = {_[0] for _ in reference_values_by_hash}
    subjects = {_[1] for _ in reference_values_by_hash}
    msg = (
        f"Plotting {len(other_runs)} runs against {run.configuration.method}"
        f" run {run.run_id} which has {len(reference_values_by_hash)} comparisons"
    )
    logger.info(msg)

    vs_count = len(other_runs)

    # Default is 6 by 6 (inches), 100dpi, making ours taller!
    # Start with a tall Figure, enough for vs_count squares plus a thin margin plot
    # Allocating a notional 5x5 for each scatter of difference plot,
    # and 1x5 or 5x1 for their histograms, and 1x5 for the spacer too
    plots_per_row = columns if columns > 0 else ceil(sqrt(vs_count))
    plots_per_col = ceil(vs_count / plots_per_row)

    done = 0
    for mode in ("scatter", "diff"):
        fig = plt.figure(figsize=(7 * plots_per_row - 1, 1 + 5 * plots_per_col))
        # Want (5+1 for plot), (1 spacer), .... (1 spacer), (5+1)
        w_ratios = tuple([5, 1] + [1, 5, 1] * (plots_per_row - 1))
        h_ratios = tuple([1] + [5] * plots_per_col)
        # Add a gridspec; adjust the subplot parameters for a square plot.
        gs = fig.add_gridspec(
            1 + plots_per_col,  # hist-x, and then one row for each comparison
            3 * plots_per_row - 1,  # plot, hist-y, spacer, ..., spacer, plot, hist-y
            width_ratios=w_ratios,
            height_ratios=h_ratios,
            left=0.15 / plots_per_row,
            right=1 - 0.15 / plots_per_row,
            bottom=0.15 / plots_per_col,
            top=1 - 0.05 / plots_per_col,
            wspace=0.05,
            hspace=0.05,
        )
        del w_ratios, h_ratios
        # Want them all to share the same x-axes, so make that first...
        scatter_axes = {
            0: fig.add_subplot(gs[1, 0]),
        }
        y_histograms = {}
        for plot_number in range(1, vs_count):
            scatter_axes[plot_number] = fig.add_subplot(
                gs[
                    1 + (plot_number // plots_per_row),
                    3 * (plot_number % plots_per_row),
                ],
                sharex=scatter_axes[0],
                sharey=scatter_axes[0] if mode == "scatter" else None,
            )
        # Add the y-histograms to the right of each plot
        for plot_number, ax in scatter_axes.items():
            y_histograms[plot_number] = fig.add_subplot(
                gs[
                    1 + (plot_number // plots_per_row),
                    1 + 3 * (plot_number % plots_per_row),
                ],
                sharey=ax,
            )
        # Adjust the appearance of all the y-value histograms
        for ax_histy in y_histograms.values():
            ax_histy.tick_params(axis="y", labelleft=False)
            ax_histy.get_xaxis().set_visible(False)
            ax_histy.spines[["top", "right", "bottom"]].set_visible(False)

        for index, ax in scatter_axes.items():
            # This doesn't catch all the bottom of a row entries
            if index // plots_per_row + 1 == plots_per_col:
                ax.set_xlabel(run.name)
            else:
                ax.tick_params(axis="x", labelbottom=False)

        # Add the top row of x-histograms
        for column in range(min(vs_count, plots_per_row)):
            ax_histx = fig.add_subplot(gs[0, column * 3], sharex=scatter_axes[0])
            ax_histx.spines[["left", "top", "right"]].set_visible(False)
            ax_histx.get_yaxis().set_visible(False)
            ax_histx.tick_params(axis="x", labelbottom=False)  # no?
            # This is a histogram of all the reference run's values - but not all will
            # have a match in any given comparison run...
            ax_histx.hist(
                reference_values_by_hash.values(),
                bins=hist_bins,
                orientation="vertical",
            )

        with Progress(*PROGRESS_BAR_COLUMNS) as progress:
            for plot_number, other_run_id in progress.track(
                enumerate(other_runs),
                description=f"Plotting {mode.ljust(8)}",
                total=len(other_runs),
            ):
                other_run = db_orm.load_run(session, other_run_id, check_complete=False)
                # Can pre-filter the query list and the suject list in the DB,
                # but ultimately need matching pairs so final filter in Python:
                other_values_by_hash = {
                    (_.query_hash, _.subject_hash): _.identity
                    for _ in other_run.comparisons()
                    .where(db_orm.Comparison.query_hash.in_(queries))
                    .where(db_orm.Comparison.subject_hash.in_(subjects))
                    .where(db_orm.Comparison.identity != None)  # noqa: E711
                    if (_.query_hash, _.subject_hash) in reference_values_by_hash
                }
                if not other_values_by_hash:
                    msg = f"Runs {run.run_id} and {other_run_id} have no comparisons in common"
                    log_sys_exit(logger, msg)
                if mode == "scatter":
                    # Don't repeat this for the diff plot
                    msg = (
                        f"Plotting {other_run.configuration.method} run {other_run_id}"
                        f" vs {run.configuration.method} run {run.run_id},"
                        f" with {len(other_values_by_hash)} comparisons in common"
                    )
                    logger.info(msg)

                # other_data dict can be smaller than ref_data!
                x_values = [
                    reference_values_by_hash[pair] for pair in other_values_by_hash
                ]
                y_values = list(other_values_by_hash.values())

                if "tsv" in formats and mode == "scatter":
                    # For simplicity (on export, but also for reuse externally),
                    # one file per pair of runs
                    filename = (
                        outdir
                        / f"{run.configuration.method}_{field}_{run.run_id}_vs_{other_run_id}.tsv"
                    )
                    with filename.open("w") as handle:
                        handle.write(f"#{run.name}\t{other_run.name}\n")
                        for x, y in zip(x_values, y_values, strict=True):
                            handle.write(f"{x}\t{y}\n")

                # Create the plot
                ax_scatter = scatter_axes[plot_number]
                ax_scatter.spines[["top", "right"]].set_visible(False)
                if mode == "diff":
                    y_values = [
                        y - x for (x, y) in zip(x_values, y_values, strict=True)
                    ]
                    # Horizontal red line y=0
                    ax_scatter.plot(
                        [min(x_values), max(x_values)], [0, 0], "-", color="r"
                    )
                else:
                    # Diagonal red y=x line
                    end_points = [
                        max(min(x_values), min(y_values)),
                        min(max(x_values), max(y_values)),
                    ]
                    ax_scatter.plot(end_points, end_points, "-", color="r")
                    del end_points

                # Draw scatter on top of red line
                ax_scatter.scatter(
                    x=x_values,
                    y=y_values,
                    s=2,
                    alpha=0.2,
                )
                ax_scatter.set_ylabel(other_run.name)
                # Now the y-value histogram
                y_histograms[plot_number].hist(
                    y_values,
                    bins=hist_bins,
                    orientation="horizontal",
                )

        for ext in formats:
            filename = (
                outdir
                / f"{run.configuration.method}_{field}_{run.run_id}_{mode}_vs_others.{ext}"
            )
            if ext == "tsv":
                pass
            else:
                fig.savefig(filename)
                done += 1
        plt.close()
    return done
