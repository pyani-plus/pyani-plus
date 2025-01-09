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

from itertools import combinations
from typing import NamedTuple

import networkx as nx  # type: ignore  # noqa: PGH003
import numpy as np  # type: ignore  # noqa: PGH003
import pandas as pd

ATTRIBUTE = "coverage"
THRESHOLD = 0.5  # default


class ClusterInfo(NamedTuple):
    """Graph structure summary."""

    n_nodes: int
    members: list
    min_cov: float | None
    min_identity: float | None
    clique: bool
    singleton: bool


def construct_complete_graph(
    cov_matrix: pd.DataFrame, id_matrix: pd.DataFrame
) -> nx.Graph:
    """Return a complete graph representing ANI results.

    Constructs an undirected complete graph for the ANI results of a given run_id.
    Nodes represent genome hashes, and edges correspond to pairwise comparisons,
    with weights assigned based on minimum coverage and average identity.
    """
    # Get list of nodes for a graph (eg. genome hashes) and a list of comparisons (eg. genome A vs genome B)
    nodes = cov_matrix.columns
    comparisons = list(combinations(nodes, 2))

    # Construct a complete graph for the results. Since the methods are not symmetrical
    # (e.g., comparisons between genome A and genome B are expected to return slightly different
    # values compared to those between genome B and genome A), we loop over the comparisons
    # and check the results in both directions. For each comparison, we use the lowest genome coverage
    # and the average genome identity to create the graph.

    rows = []
    for genome1, genome2 in comparisons:
        datadict = {
            "genome1": genome1,
            "genome2": genome2,
            "coverage": min(cov_matrix[genome1][genome2], cov_matrix[genome2][genome1]),
            "identity": np.mean(
                [id_matrix[genome1][genome2], id_matrix[genome2][genome1]]
            ),
        }
        rows.append(datadict)

    node_data = pd.DataFrame(rows)

    return nx.from_pandas_edgelist(
        node_data, "genome1", "genome2", ["coverage", "identity"]
    )


def remove_edges_below_threshold(
    graph: nx.Graph, attribute: str | None, threshold: float | None
) -> nx.Graph:
    """Return a graph where edges with weights of the specified attribute fall below the given threshold."""
    edges_to_remove = [
        (node1, node2)
        for node1, node2, attrs in graph.edges(data=True)
        if attrs[attribute] < threshold
    ]
    graph.remove_edges_from(edges_to_remove)

    return graph


def analyse_subgraphs(graph: nx.Graph) -> list[ClusterInfo]:
    """Analyse subgraphs in a graph and classify them into cliques, singletons, or non cliques."""
    new_graph = graph.copy()
    subgraphs = [graph.subgraph(comp) for comp in nx.connected_components(new_graph)]
    node_degrees = dict(graph.degree)

    subgraphs_info = []
    for subgraph in subgraphs:
        nodes = sorted(subgraph.nodes)
        n_nodes = len(nodes)
        edges = list(subgraph.edges(data=True))

        # Check if all nodes are fully connected (clique)
        is_clique = all(node_degrees[node] == n_nodes - 1 for node in nodes)
        # Check subgraph is a siglenton
        is_singleton = n_nodes == 1

        # Get minimum attributes (eg. coverage and identity) if edges exist
        min_cov = min(
            (attrs["coverage"] for node1, node2, attrs in edges), default=None
        )
        min_identity = min(
            (attrs["identity"] for node1, node2, attrs in edges), default=None
        )

        subgraphs_info.append(
            ClusterInfo(
                n_nodes=n_nodes,
                members=nodes,
                min_cov=min_cov,
                min_identity=min_identity,
                clique=is_clique and not is_singleton,
                singleton=is_singleton,
            )
        )

    return subgraphs_info
