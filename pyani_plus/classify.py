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

from collections.abc import Callable
from itertools import combinations
from pathlib import Path
from typing import NamedTuple

import networkx as nx
import numpy as np
import pandas as pd
from rich.progress import Progress

from pyani_plus import PROGRESS_BAR_COLUMNS

AGG_FUNCS = {
    "min": min,
    "max": max,
    "mean": np.mean,
}

MIN_COVERAGE = 0.50


class CliqueInfo(NamedTuple):
    """Graph structure summary."""

    n_nodes: int
    min_cov: float | None
    min_identity: float | None
    members: list


def construct_graph(
    cov_matrix: pd.DataFrame,
    id_matrix: pd.DataFrame,
    coverage_agg: Callable,
    identity_agg: Callable,
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
        identity = identity_agg(
            [id_matrix[genome1][genome2], id_matrix[genome2][genome1]]
        )
        # Add edge only if both coverage and identity are valid
        if pd.notna(coverage) and pd.notna(identity) and coverage > min_coverage:
            graph.add_edge(genome1, genome2, coverage=coverage, identity=identity)

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
    for component in connected_components:
        subgraph = graph.subgraph(component).copy()
        if is_clique(subgraph):  # Check if the subgraph is a clique
            cliques.append(subgraph)

    return cliques


def find_cliques_recursively(
    graph: nx.Graph,
    progress=None,  # noqa: ANN001
    task=None,  # noqa: ANN001
) -> list:
    """Return all unique cliques within a set of genomes based on ANI results.

    These cliques are identified recursively by iteratively removing edges
    with the lowest identity at each step and analysing the resulting subgraphs.
    """
    cliques: list(nx.Graph) = []  # type: ignore  # noqa: PGH003

    # If the graph has only one node stop recursion and return list of cliques
    if len(graph.nodes) == 1:
        cliques.append(graph)
        return cliques
    # Recording cliques
    if is_clique(graph):
        cliques.append(graph.copy())

    edges = graph.edges(data=True)
    edges = sorted(edges, key=lambda edge: edge[2]["identity"])

    # Initialise the progress bar only at the top level
    if progress is None:
        with Progress(*PROGRESS_BAR_COLUMNS) as progress:  # noqa: PLR1704
            task = progress.add_task("Processing edges...", total=len(edges))
            return find_cliques_recursively(graph, progress, task)

    # Remove edges with the lowest weight and identify cliques
    while edges:
        edge_to_remove = edges.pop(0)
        progress.update(task, advance=1)  # Update the progress bar
        graph.remove_edge(edge_to_remove[0], edge_to_remove[1])

        connected_components = list(nx.connected_components(graph))
        if len(connected_components) > 1:
            for component in connected_components:
                subgraph = graph.subgraph(component).copy()
                cliques.extend(find_cliques_recursively(subgraph, progress, task))
            return cliques

    return cliques


def compute_classify_output(
    cliques: list, method: str, outdir: Path
) -> list[CliqueInfo]:
    """Return list of CliqueInfo describing all cliques found and save them to .tsv file."""
    clique_data = [
        CliqueInfo(
            n_nodes=len(clique.nodes),
            min_cov=min(
                (attrs["coverage"] for _, _, attrs in clique.edges(data=True)),
                default=None,
            ),
            min_identity=min(
                (attrs["identity"] for _, _, attrs in clique.edges(data=True)),
                default=None,
            ),
            members=list(clique.nodes),
        )
        for clique in cliques
    ]

    clique_df = pd.DataFrame(clique_data)
    clique_df["members"] = clique_df["members"].apply(lambda x: ",".join(x))

    output_file = outdir / f"{method}_classify.tsv"
    # Round coverage and identity values to 7 decimal places before saving
    clique_df.round(7).to_csv(output_file, sep="\t", index=False)

    return clique_data
