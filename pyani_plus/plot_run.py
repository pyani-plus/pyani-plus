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

import warnings

import seaborn as sns
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
from pandas import DataFrame

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


def plot_heatmap(
    matrix: DataFrame,
    output_stem: str,
    color_scheme: str,
    formats: tuple[str, ...] = ("tsv", "png", "jpg", "svg", "pdf"),
) -> None:
    """Plot the given matrix as a heatmap."""
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
        filename = output_stem + "." + ext
        if ext == "tsv":
            # Apply the clustering reordering to match the figure:
            matrix = matrix.iloc[
                figure.dendrogram_row.reordered_ind,
                figure.dendrogram_row.reordered_ind,
            ]
            matrix.to_csv(filename, sep="\t")
        else:
            figure.savefig(filename)
