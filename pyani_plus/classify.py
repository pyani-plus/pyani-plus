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
"""Code to implement the classify method indetnded to identify cliques within a set of genomes."""

import math
from collections import defaultdict
from collections.abc import Callable
from itertools import combinations
from pathlib import Path
from typing import NamedTuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import cm, patches
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from rich.progress import Progress

from pyani_plus import GRAPHICS_FORMATS, PROGRESS_BAR_COLUMNS
from pyani_plus.public_cli_args import EnumModeClassify

AGG_FUNCS = {
    "min": min,
    "max": max,
    "mean": np.mean,
}

MIN_COVERAGE = 0.50

MODE = EnumModeClassify.identity  # constant for CLI default


class CliqueInfo(NamedTuple):
    """Graph structure summary."""

    n_nodes: int
    max_cov: float | None
    min_score: float | None
    max_score: float | None
    members: list


def construct_graph(
    cov_matrix: pd.DataFrame,
    score_matrix: pd.DataFrame,
    coverage_agg: Callable,
    score_agg: Callable,
    min_coverage: float,
) -> nx.Graph:
    """Return a graph representing ANI results.

    Constructs an undirected graph for the ANI results of a given run_id.
    Nodes represent genome, and edges correspond to pairwise comparisons,
    with weights assigned based on the minimum coverage and the average identity
    (by default). Edges are not added between genome pairs if no alignment is
    found or if the coverage is below 50%.
    """
    # Create an empty graph
    graph = nx.Graph()

    # Get list of nodes for a graph (eg. genome hashes) and a list of comparisons (eg. genome A vs genome B)
    nodes = cov_matrix.columns
    graph.add_nodes_from(nodes)
    comparisons = list(combinations(nodes, 2))

    # Add edges to the graph based on ANI results. Since the methods are not symmetrical
    # (e.g., comparisons between genome A and genome B are expected to return slightly different
    # values compared to those between genome B and genome A), we loop over the comparisons
    # and check the results in both directions. For each comparison, by default, we use
    # the lowest genome coverage and the average genome identity to create the graph.
    # However, we allow the end-user to choose alternatives, such as min or max. No edges
    # are added if no alignment is found.
    for genome1, genome2 in comparisons:
        coverage = coverage_agg(
            [cov_matrix[genome1][genome2], cov_matrix[genome2][genome1]]
        )
        score = score_agg(
            [score_matrix[genome1][genome2], score_matrix[genome2][genome1]]
        )
        # Add edge only if both coverage and identity are valid
        if pd.notna(coverage) and pd.notna(score) and coverage > min_coverage:
            graph.add_edge(genome1, genome2, coverage=coverage, score=score)

    return graph


def is_clique(graph: nx.Graph) -> bool:
    """Return True if the subgraph is a clique."""
    n_nodes = len(graph.nodes)
    return len(graph.edges) == n_nodes * (n_nodes - 1) / 2


def find_initial_cliques(graph: nx.Graph) -> list:
    """Return all unique cliques in the given graph.

    Since the initial graph has edges removed below the 50% coverage
    threshold (by default), it is possible that the graph contains
    subgraphs that are potential cliques, which we may want to identify
    before removing any edges.
    """
    cliques: list(nx.Graph) = []  # type: ignore  # noqa: PGH003
    connected_components = list(nx.connected_components(graph))
    edges = nx.get_edge_attributes(graph, "score")

    identity = min(edges.values()) if edges else None
    for component in connected_components:
        subgraph = graph.subgraph(component).copy()
        if is_clique(subgraph):  # Check if the subgraph is a clique
            cliques.append((subgraph, identity))

    return cliques


def find_cliques_recursively(
    graph: nx.Graph,
    progress=None,  # noqa: ANN001
    task=None,  # noqa: ANN001
    min_score: float | None = None,
) -> list[tuple]:
    """Return all cliques within a set of genomes based on ANI
    results, along with the identity edge that formed each clique, as a tuple.

    These cliques are identified recursively by iteratively removing edges
    with the lowest identity at each step and analysing the resulting subgraphs.
    """  # noqa: D205
    cliques = []

    # If the graph has only one node stop recursion and return list of cliques
    if len(graph.nodes) == 1:
        cliques.append((graph, min_score))
        return cliques
    # Recording cliques
    if is_clique(graph):
        cliques.append((graph.copy(), min_score))

    edges = graph.edges(data=True)
    edges = sorted(edges, key=lambda edge: edge[2]["score"])

    # Initialise the progress bar only at the top level
    if progress is None:
        with Progress(*PROGRESS_BAR_COLUMNS) as progress:  # noqa: PLR1704
            task = progress.add_task("Processing edges...", total=len(edges))
            return find_cliques_recursively(graph, progress, task, min_score=min_score)

    # Remove edges with the lowest weight, identify cliques and retain the identity edge that formed each clique.

    while edges:
        edge_to_remove = edges.pop(0)
        break_edge = (edge_to_remove[0], edge_to_remove[1])
        min_score = graph.get_edge_data(*break_edge).get("score")
        progress.update(task, advance=1)  # Update the progress bar
        graph.remove_edge(edge_to_remove[0], edge_to_remove[1])

        connected_components = list(nx.connected_components(graph))
        if len(connected_components) > 1:
            for component in connected_components:
                subgraph = graph.subgraph(component).copy()
                cliques.extend(
                    find_cliques_recursively(
                        subgraph,
                        progress,
                        task,
                        min_score=min_score,
                    )
                )
            return cliques

    return cliques


def get_unique_cliques(
    initial_cliques: list[tuple], recursive_cliques: list[tuple]
) -> list[tuple]:
    """Return only unique cliques, along with the identity edge that formed each clique, as a tuple."""
    unique_cliques = {
        frozenset(graph.nodes): (graph, edge) for graph, edge in initial_cliques
    }
    unique_cliques.update(
        {
            frozenset(graph.nodes): (graph, edge)
            for graph, edge in recursive_cliques
            if frozenset(graph.nodes) not in unique_cliques
        }
    )

    return list(unique_cliques.values())


def get_genome_cligue_ids(dataframe: pd.DataFrame, suffix: str) -> dict:
    """Return a dictionary mapping each genome to a list of clique ids it belongs to."""
    # Members of each clique are split into separate entries.
    # None in max_{score} indicated a singleton which will be 1.0 for identity and 0.0 for tANI
    dataframe[f"max_{suffix}"] = dataframe[f"max_{suffix}"].fillna(
        1.0 if suffix == "identity" else 0.0
    )
    dataframe["members"] = dataframe["members"].str.split(",")

    genome_clique_ids = defaultdict(list)
    for idx, row in dataframe.iterrows():
        for genome in row["members"]:
            genome_clique_ids[genome].append(idx)

    return genome_clique_ids


def get_genome_order(genome_clique_ids: dict) -> dict:
    """Return dictionary mapping each genome to a position on the y-axis."""
    # Sort genome_clique_ids by the clique IDs and map each genome to a unique y-position.
    sorted_genomes = sorted(
        genome_clique_ids, key=lambda genome: genome_clique_ids[genome]
    )
    return {genome: idx for idx, genome in enumerate(sorted_genomes)}


def plot_classify(  # noqa: PLR0913, PLR0915
    genome_positions: dict,
    dataframe: pd.DataFrame,
    outdir: Path,
    method: str,
    score: str,
    vertical_line: float,
    formats: tuple[str, ...] = GRAPHICS_FORMATS,
) -> None:
    """Plot the classify results for a given run.

    The function generates a plot with 4 vertically stacked subplots:

    1. ax1: The number of genomes in cliques and as singletons at different %identity.
    2. ax2: The percentage of all genomes that are found in cliques and as singletons at different %identity.
    3. ax3: The lifespan (range of identity) of each clique.
    4. ax4: Colorbar/Legend for the colours of cliques in ax3.
    """
    # Set up initial plot (eg. figsize, subplots and colormap)
    num_genomes = len(genome_positions)

    # Dynamically adjust figure height based on the number of genomes
    label_spacing = 0.15
    min_height = 15  # Prevent plot from being too small
    height = max(num_genomes * label_spacing, min_height)

    # Dynamically adjust adjust y-axis label font size
    font_size = max(6, min(12, 300 // num_genomes))

    # Limit hspace to a maximum of 0.1
    hspace = min(0.1, 10 / num_genomes)
    # Create a figure with 3 vertically stacked subplots
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(
        4,
        1,
        figsize=(15, height),
        gridspec_kw={
            "height_ratios": [0.7, 0.7, max(5, num_genomes * 0.1), 0.2],
            "hspace": hspace,
        },
        sharex=True,
    )
    fig.subplots_adjust(left=0.2, top=0.85, hspace=hspace)

    # Ensure x-axis tick labels on all plots
    ax1.tick_params(labelbottom=True)
    ax2.tick_params(labelbottom=True)

    # Setting a colorbar on the bottom of the plot for cliques based on min_identity
    norm = Normalize(
        vmin=math.floor(dataframe[f"min_{score}"].min() * 100) / 100 - 0.01,
        vmax=dataframe[f"min_{score}"].max(),
    )
    # Create a colormap with alpha=0.7
    cmap_hot = cm.hot

    # ----- The number of genomes in cliques and as singletons at different %identity (ax1) -----
    # Define bins for identity range
    identity_bins = np.linspace(
        math.floor(dataframe[f"min_{score}"].min() * 100) / 100, 1.0, 100
    )
    clique_counts = np.zeros_like(
        identity_bins[:-1]
    )  # Stores counts of genomes in cliques
    singleton_counts = np.zeros_like(
        identity_bins[:-1]
    )  # Stores counts of singleton genomes

    # Iterate through each row in the dataframe to populate count arrays
    for _, row in dataframe.iterrows():
        x_start, x_end = row[f"min_{score}"], row[f"max_{score}"]  # Identity range
        mask = (identity_bins[:-1] >= x_start) & (
            identity_bins[:-1] <= x_end
        )  # Select bins in range

        if len(row["members"]) > 1:
            clique_counts[mask] += len(row["members"])  # Increment clique genome count
        else:
            singleton_counts[mask] += 1  # Increment singleton count

    # Plot the number of genomes in cliques
    ax1.plot(
        identity_bins[:-1],
        clique_counts,
        color="blue",
        linewidth=2,
        label="Genomes in Cliques",
    )
    ax1.fill_between(identity_bins[:-1], clique_counts, color="blue", alpha=0.3)

    # Plot the number of singleton genomes
    ax1.plot(
        identity_bins[:-1],
        singleton_counts,
        color="red",
        linewidth=2,
        linestyle="--",
        label="Singleton Genomes",
    )
    # Customise ax1: set y-label, grid style, and display legend
    ax1.set_ylabel("Number of \n Genomes", fontsize=10)
    ax1.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)  # noqa: FBT003
    ax1.legend()

    # --- The percentage of all genomes that are found in cliques and as singletons at different %identity (ax2) ---
    total_genomes = clique_counts + singleton_counts  # Total count per identity bin
    percent_total_genomes = (total_genomes / num_genomes) * 100  # Convert to percentage

    # Plot percentage of genomes (cliques + singletons)
    ax2.plot(
        identity_bins[:-1],
        percent_total_genomes,
        color="green",
        linewidth=2,
        label="% Genomes",
    )
    ax2.fill_between(
        identity_bins[:-1], percent_total_genomes, color="green", alpha=0.3
    )
    # Customise ax2: set y-label, grid style, and display legend
    ax2.set_ylabel("Percentage of \n All Genomes", fontsize=10)
    ax2.set_ylim(0, 100)
    ax2.grid(True, linestyle="--", linewidth=0.5, alpha=0.9)  # noqa: FBT003
    ax2.legend()

    # --- The lifespan (range of identity) of each clique (ax3) ---
    for _, row in dataframe.iterrows():
        x_start, x_end = row[f"min_{score}"], row[f"max_{score}"]
        y_positions = [genome_positions[g] for g in row["members"]]
        y_min, y_max = min(y_positions), max(y_positions)

        if len(row["members"]) == 1:
            # Plot single-member groups as dashed horizontal lines
            ax3.hlines(
                y=y_min,
                xmin=x_start,
                xmax=x_end,
                colors="grey",
                linestyles="dashed",
                linewidth=1.5,
            )
        else:
            # Plot cliques as rectangles with heatmap coloring
            color = cmap_hot(norm(x_start))
            rect = patches.Rectangle(
                (x_start, y_min - 0.4),
                x_end - x_start,
                y_max - y_min + 0.8,
                linewidth=1,
                edgecolor="black",
                facecolor=color,
                alpha=0.8,
            )
            ax3.add_patch(rect)

    ax3.set_xlabel(f"{score}")
    ax3.set_ylabel("Genomes", fontsize=6)
    ax3.set_yticks(range(num_genomes))
    ax3.set_yticklabels(genome_positions.keys(), fontsize=font_size)
    ax3.yaxis.set_label_position("right")
    ax3.yaxis.tick_right()
    ax3.set_xlim(
        math.floor(dataframe[f"min_{score}"].min() * 100) / 100 - 0.01,
        dataframe[f"max_{score}"].max(),
    )
    ax3.set_ylim(-1, num_genomes)
    x_value = (
        vertical_line
        if vertical_line != 0.95  # noqa: PLR2004
        else (0.95 if score == "identity" else -0.323)
    )
    ax3.axvline(x=x_value, color="red", linewidth=2, linestyle="--")
    ax3.grid(True, linestyle="--", linewidth=0.5, alpha=0.9)  # noqa: FBT003

    # -------- Colorbar (ax4) ----------
    gradient_values = np.linspace(norm.vmin, norm.vmax, 10000)
    lines = [[(value, 0), (value, 1)] for value in gradient_values]
    colors = cmap_hot(norm(gradient_values))

    line_collection = LineCollection(lines, colors=colors, linewidths=0.5)
    ax4.add_collection(line_collection)
    ax4.set_xlim(norm.vmin, norm.vmax)
    ax4.set_ylim(0, 1)
    ax4.set_xlabel(f"Min {score}", fontsize=10)
    ax4.xaxis.set_label_position("bottom")
    ax4.set_yticks([])  # No y-axis
    ax4.tick_params(axis="x", labelsize=10, direction="out")
    # Save the figure in all standard graphics formats
    for ext in formats:
        if ext != "tsv":
            plt.savefig(
                outdir / f"{method}_classify_plot.{ext}",
                format=ext,
                bbox_inches="tight",
            )


def compute_classify_output(
    cliques: list, method: str, outdir: Path, column_map: dict
) -> tuple[list[CliqueInfo], pd.DataFrame]:
    """Return list of CliqueInfo describing all cliques found and save them to .tsv file."""
    clique_data = [
        CliqueInfo(
            n_nodes=len(clique.nodes),
            max_cov=min(
                (attrs["coverage"] for _, _, attrs in clique.edges(data=True)),
                default=None,
            ),
            min_score=edge_form,
            max_score=min(
                (attrs["score"] for _, _, attrs in clique.edges(data=True)),
                default=None,
            ),
            members=list(clique.nodes),
        )
        for clique, edge_form in cliques
    ]

    clique_df = pd.DataFrame(clique_data)
    clique_df["members"] = clique_df["members"].apply(lambda x: ",".join(x))

    # Rename columns based on mode
    clique_df = clique_df.rename(columns=column_map)

    output_file = outdir / f"{method}_classify.tsv"
    # Round coverage and identity values to 7 decimal places before saving
    clique_df.round(7).to_csv(output_file, sep="\t", index=False)

    return clique_data, clique_df
