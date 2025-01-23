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
"""Tests for the classify implementation.

These tests are intended to be run from the repository root using:

pytest -v
"""

import networkx as nx  # type: ignore  # noqa: PGH003
import networkx.algorithms.isomorphism as iso
import numpy as np
import pandas as pd
import pytest

from pyani_plus import classify


@pytest.fixture
def known_graph_with_dataframes() -> tuple[nx.Graph, pd.DataFrame, pd.DataFrame]:
    """Return nx.Graph with known weights (coverage and identity) and corresponding coverage and identity DataFrames.

    The graph:
    - Contains 5 nodes named 'genome_1' to 'genome_5'.
    - Each pair of nodes is connected by an edge with the following attributes:
        - 'coverage': Fixed at 1.0 for all edges.
        - 'identity': Starts at 0.80 and increases by 0.01 for each subsequent edge.
    - Edges are added in a deterministic order such that the identity values remain consistent
      across multiple runs.
    """
    # Create an empty graph and add nodes
    graph = nx.Graph()
    nodes = [f"genome_{i}" for i in range(1, 6)]
    graph.add_nodes_from(nodes)

    # Create empty DataFrames
    coverage_df = pd.DataFrame(None, index=nodes, columns=nodes)
    identity_df = pd.DataFrame(None, index=nodes, columns=nodes)

    # Add edges with attributes
    identity_value = 0.80
    edges = []
    for i, node1 in enumerate(nodes):
        for node2 in nodes[i + 1 :]:
            # Add edges to the graph
            edge_attributes = {"coverage": 1.0, "identity": round(identity_value, 2)}
            edges.append((node1, node2, edge_attributes))

            # Update DataFrames ()
            coverage_df.loc[node1, node2] = 1.0
            coverage_df.loc[node2, node1] = 1.0
            identity_df.loc[node1, node2] = round(identity_value, 2)
            identity_df.loc[node2, node1] = round(identity_value, 2)

            # Increment identity value
            identity_value += 0.01

    graph.add_edges_from(edges)

    return graph, coverage_df, identity_df


@pytest.fixture
def expected_cliques(
    known_graph_with_dataframes: tuple[nx.Graph, pd.DataFrame, pd.DataFrame],
) -> list:
    """Generate a list of expected cliques based on the known graph.

    Since, the known graph is constructed such that nodes are connected in order,
    and edge 'identity' values increment by 0.01, starting at 0.80. Removing
    'genome_1' first results in the removal of edges with the lowest identity
    values, leaving behind the first expected clique. We can generate a list
    of all possible cliques that should be identified by the classification method
    """
    # Copy the original graph to avoid modifying the fixture
    original_graph = known_graph_with_dataframes[0].copy()
    graph = original_graph.copy()

    # Define empty list to which we will append known cliques
    known_cliques = []

    # Iterate through nodes in the graph
    for node in list(graph.nodes):
        if node != "genome_5":
            # Add single-node graph as a clique
            single_node_clique = nx.Graph()
            single_node_clique.add_node(node)
            known_cliques.append(single_node_clique)

            # Remove the node and add the resulting graph as a clique
            graph.remove_node(node)
            known_cliques.append(graph.copy())  # Copy to preserve graph state

    # Include the initial graph as the first clique
    return [original_graph, *known_cliques]


def test_construct_graph(
    known_graph_with_dataframes: tuple[nx.Graph, pd.DataFrame, pd.DataFrame],
) -> None:
    """Check construction of the initial graph."""
    graph, coverage_df, identity_df = known_graph_with_dataframes

    # Comparison function for a numerical edge attribute.
    edge_match = iso.numerical_edge_match("coverage", "identity")

    # Check the isomorphism of a graph with the edge_match function
    assert nx.is_isomorphic(
        graph,
        classify.construct_graph(coverage_df, identity_df, min, np.mean, 0.5),
        edge_match=edge_match,
    )


def test_is_clique(
    known_graph_with_dataframes: tuple[nx.Graph, pd.DataFrame, pd.DataFrame],
) -> None:
    """Check cliques are identified."""
    graph = known_graph_with_dataframes[0]
    assert classify.is_clique(graph) is True


def test_find_initial_cliques(
    known_graph_with_dataframes: tuple[nx.Graph, pd.DataFrame, pd.DataFrame],
) -> None:
    """Check all possible cliques are identified in the initial iteration."""
    graph = known_graph_with_dataframes[0]
    found_cliques = classify.find_initial_cliques(graph)
    edge_match = iso.numerical_edge_match("coverage", "identity")

    # Get the connected components as subgraphs
    connected_components = [
        graph.subgraph(component).copy() for component in nx.connected_components(graph)
    ]

    # Check the number of identified cliques
    assert len(connected_components) == len(found_cliques), "Clique count mismatch"

    # Check the clique structures
    for expected_clique, found_clique in zip(
        connected_components, found_cliques, strict=False
    ):
        assert nx.is_isomorphic(expected_clique, found_clique, edge_match=edge_match), (
            "Clique structure mismatch"
        )


def test_find_cliques(
    known_graph_with_dataframes: tuple[nx.Graph, pd.DataFrame, pd.DataFrame],
    expected_cliques: list,
) -> None:
    """Check all possible cliques are identified in the initial iteration."""
    graph = known_graph_with_dataframes[0]
    found_cliques = classify.find_cliques_recursively(graph)
    edge_match = iso.numerical_edge_match("coverage", "identity")

    # Check the number of identified cliques
    assert len(found_cliques) == len(expected_cliques), "Clique count mismatch"

    for expected_clique, found_clique in zip(
        expected_cliques, found_cliques, strict=False
    ):
        assert nx.is_isomorphic(expected_clique, found_clique, edge_match=edge_match), (
            "Clique structure mismatch"
        )
