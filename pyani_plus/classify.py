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

AGG_FUNCS = {
    "min": min,
    "max": max,
    "mean": np.mean,
}


class CliqueInfo(NamedTuple):
    """Graph structure summary."""

    members: list
    n_nodes: int
    min_cov: float | None
    min_identity: float | None


def construct_graph(
    cov_matrix: pd.DataFrame,
    id_matrix: pd.DataFrame,
    coverage_agg: Callable,
    identity_agg: Callable,
) -> nx.Graph:
    """Return a graph representing ANI results.

    Constructs an undirected graph for the ANI results of a given run_id.
    Nodes represent genome, and edges correspond to pairwise comparisons,
    with weights assigned based on the minimum coverage and the average identity
    (by default). No edges are added between pairs of genomes when no alignment
    is found.
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
        if pd.notna(coverage) and pd.notna(identity):
            graph.add_edge(genome1, genome2, coverage=coverage, identity=identity)

    return graph


def is_clique(graph: nx.Graph) -> bool:
    """Return True if the subgraph is a clique."""
    nodes = list(graph.nodes)
    return all(graph.degree[node] == len(nodes) - 1 for node in nodes)


def remove_lowest_edge(subgraph: nx.Graph, attribute: str) -> nx.Graph:
    """Remove the lowest edge from the subgraph."""
    lowest_edge = min([attrs[attribute] for n1, n2, attrs in subgraph.edges(data=True)])

    edges_to_remove = [
        (n1, n2)
        for n1, n2, attrs in subgraph.edges(data=True)
        if attrs[attribute] == lowest_edge
    ]
    subgraph.remove_edges_from(edges_to_remove)

    return subgraph


def find_cliques(
    graph: nx.Graph, attribute: str, seen_cliques: list | None = None
) -> list:
    """Return all unique cliques within a set of genomes based on ANI results."""
    if seen_cliques is None:
        seen_cliques = []

    cliques: list(nx.Graph) = []  # type: ignore  # noqa: PGH003

    # If the graph has only one node stop recursion and return list of cliques
    if len(graph.nodes) == 1:
        return cliques

    # This is up for discussion, but at the moment, the code will return only unique cliques and
    # will not report the same cliques twice. Even if they appear multiple times across different
    # levels of recursion. This is achieved by only appending cliques that have not been
    # previously seen in `seen_cliques`. NOTE: We might want to keep track of all cliques to monitor
    # how long a specific clique stays intact (e.g., which range of thresholds the clique remains
    # stable at). But, this can also be done by examining the edge attributes/weights.
    if is_clique(graph):
        clique_repr = (tuple(sorted(graph.nodes)), tuple(sorted(graph.edges)))
        if clique_repr not in seen_cliques:
            cliques.append(graph.copy())
            seen_cliques.append(clique_repr)

    temp_graph = graph.copy()

    # Remove edges one by one and recursively process connected components
    while temp_graph.number_of_edges() > 0:
        remove_lowest_edge(temp_graph, attribute)

        # Check for cliques in connected components
        for subgraph_nodes in nx.connected_components(temp_graph):
            subgraph = temp_graph.subgraph(subgraph_nodes).copy()
            cliques.extend(find_cliques(subgraph, attribute, seen_cliques))

    return cliques


def compute_classify_output(
    cliques: list, method: str, attribute: str, outdir: Path
) -> list[CliqueInfo]:
    """Return list of CliqueInfo describing all cliques found and save them to .tsv file."""
    clique_data = [
        CliqueInfo(
            members=list(clique.nodes),
            n_nodes=len(clique.nodes),
            min_cov=min(
                (attrs["coverage"] for _, _, attrs in clique.edges(data=True)),
                default=None,
            ),
            min_identity=min(
                (attrs["identity"] for _, _, attrs in clique.edges(data=True)),
                default=None,
            ),
        )
        for clique in cliques
    ]

    clique_df = pd.DataFrame(clique_data)
    clique_df["members"] = clique_df["members"].apply(lambda x: ",".join(x))

    output_file = outdir / f"{method}_{attribute}_classify.tsv"
    clique_df.to_csv(output_file, sep="\t", index=False)

    return clique_data
