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
def one_node_no_edges_graph() -> nx.Graph:
    """Return graph with one node and no edges and the expected identity edge that formed the clique."""
    graph = nx.Graph()
    graph.add_node("genome_1")
    return (graph, None)


@pytest.fixture
def two_nodes_no_edges_graph() -> nx.Graph:
    """Return graph with two nodes and no edges."""
    graph = nx.Graph()
    graph.add_nodes_from(["genome_1", "genome_2"])
    return graph


@pytest.fixture
def two_nodes_no_edges_cliques() -> list[nx.Graph]:
    """Return list of all possible cliques for initial graph with two nodes and no edges
    and the expected identity edge that formed the clique.
    """  # noqa: D205
    clique_1 = nx.Graph()
    clique_1.add_node("genome_1")

    clique_2 = nx.Graph()
    clique_2.add_node("genome_2")

    return [(clique_1, None), (clique_2, None)]


@pytest.fixture
def two_nodes_no_edges_dataframes() -> list[pd.DataFrame]:
    """Return dataframes for building initial graph with two nodes and no edges."""
    # Define the genomes
    genomes = ["genome_1", "genome_2"]

    # Create the DataFrame directly with initial values
    # When building graph with `construct_graph` no edges should be added if coverage is <0.50
    coverage_df = pd.DataFrame(
        [[1.0, 0.40], [0.40, 1.0]], index=genomes, columns=genomes
    )
    identity_df = pd.DataFrame(
        [[1.0, 0.80], [0.80, 1.0]], index=genomes, columns=genomes
    )

    return [coverage_df, identity_df]


@pytest.fixture
def two_nodes_one_edge_graph() -> nx.Graph:
    """Return graph with two nodes and one edge and the expected identity edge that formed the clique."""
    graph = nx.Graph()
    graph.add_edge("genome_1", "genome_2", score=0.999310, coverage=0.6774176803)
    return (graph, 0.999310)


@pytest.fixture
def two_nodes_one_edge_cliques() -> list[tuple[nx.Graph, float | str]]:
    """Return list of all possible cliques with initial graph of two nodes and one edge
    and the expected identity edge that formed each clique.
    """  # noqa: D205
    clique_1 = nx.Graph()
    clique_1.add_edge("genome_1", "genome_2", score=0.999310, coverage=0.6774176803)

    clique_2 = nx.Graph()
    clique_2.add_node("genome_1")

    clique_3 = nx.Graph()
    clique_3.add_node("genome_2")

    return [(clique_1, 0.999310), (clique_2, 0.999310), (clique_3, 0.999310)]


@pytest.fixture
def known_complex_graph() -> nx.Graph:
    """Return nx.Graph with six nodes and known weights (coverage and identity).

    This six node graph first splits into two cliques of three-nodes:
    - genome_1, genome_5 and genome_6
    - genome_2, genome_3 and genome 4

    Each three-node clique then splits into a two-node clique:
    - genome_2 and genome_3
    - genome_1 and genome_6

    Finally, all nodes become singletons.
    """
    graph = nx.Graph()

    graph.add_edge("genome_1", "genome_2", score=0.85, coverage=1.0)
    graph.add_edge("genome_1", "genome_5", score=0.96, coverage=1.0)
    graph.add_edge("genome_1", "genome_6", score=0.99, coverage=1.0)
    graph.add_edge("genome_2", "genome_3", score=0.97, coverage=1.0)
    graph.add_edge("genome_2", "genome_4", score=0.967, coverage=1.0)
    graph.add_edge("genome_3", "genome_4", score=0.95, coverage=1.0)
    graph.add_edge("genome_4", "genome_5", score=0.86, coverage=1.0)
    graph.add_edge("genome_5", "genome_6", score=0.98, coverage=1.0)

    return graph


@pytest.fixture
def known_complex_cliques() -> list[nx.Graph]:
    """Return a list of nx.Graph objects representing known cliques with coverage and
    identity attributes and the expected identity edge that formed each clique.
    """  # noqa: D205
    # Define cliques
    cliques_data: list[list] = [
        # Clique 1
        [
            ("genome_1", "genome_5", {"score": 0.96, "coverage": 1.0}),
            ("genome_1", "genome_6", {"score": 0.99, "coverage": 1.0}),
            ("genome_5", "genome_6", {"score": 0.98, "coverage": 1.0}),
        ],
        # Clique 2
        [("genome_1", "genome_6", {"score": 0.99, "coverage": 1.0})],
        # Clique 3
        [("genome_1",)],
        # Clique 4
        [("genome_6",)],
        # Clique 5
        [("genome_5",)],
        # Clique 6
        [
            ("genome_2", "genome_3", {"score": 0.97, "coverage": 1.0}),
            ("genome_2", "genome_4", {"score": 0.967, "coverage": 1.0}),
            ("genome_3", "genome_4", {"score": 0.95, "coverage": 1.0}),
        ],
        # Clique 7
        [("genome_2", "genome_3", {"score": 0.97, "coverage": 1.0})],
        # Clique 8
        [("genome_2",)],
        # Clique 9
        [("genome_3",)],
        # Clique 10
        [("genome_4",)],
    ]

    edges = [0.86, 0.98, 0.99, 0.99, 0.98, 0.86, 0.967, 0.97, 0.97, 0.967]

    # Create graph objects for each clique
    cliques = []
    for clique_data, edge_value in zip(cliques_data, edges, strict=False):
        graph = nx.Graph()
        for edge in clique_data:
            if len(edge) == 1:  # Single node
                graph.add_node(edge[0])
            else:  # Edge with attributes
                graph.add_edge(edge[0], edge[1], **edge[2])
        cliques.append((graph, edge_value))

    return cliques


def test_construct_graph(
    two_nodes_no_edges_graph: nx.Graph,
    two_nodes_no_edges_dataframes: list[pd.DataFrame],
) -> None:
    """Check construction of the initial graph."""
    coverage_df, identity_df = two_nodes_no_edges_dataframes

    # Comparison function for a numerical edge attribute.
    edge_match = iso.numerical_edge_match("coverage", "identity")

    # Check the isomorphism of a graph with the edge_match function
    assert nx.is_isomorphic(
        two_nodes_no_edges_graph,
        classify.construct_graph(coverage_df, identity_df, min, np.mean, 0.5),
        edge_match=edge_match,
    )


def test_is_clique(
    two_nodes_one_edge_graph: nx.Graph,
) -> None:
    """Check cliques are identified."""
    assert classify.is_clique(two_nodes_one_edge_graph[0]) is True


def test_find_initial_cliques(two_nodes_one_edge_graph: tuple[nx.Graph, str]) -> None:
    """Check all possible cliques and their associated values are identified in the initial iteration."""
    graph, edge_value = two_nodes_one_edge_graph
    found_cliques = classify.find_initial_cliques(graph)
    edge_match = iso.numerical_edge_match("coverage", "score")

    # Get the connected components as subgraphs
    connected_components = [
        (graph.subgraph(component).copy(), edge_value)
        for component in nx.connected_components(graph)
    ]

    # Check the number of identified cliques
    assert len(connected_components) == len(found_cliques), "Clique count mismatch"

    # Check that both the graph structures and identity break edges values match
    for (expected_clique, expected_edge), (found_clique, found_edge) in zip(
        connected_components, found_cliques, strict=False
    ):
        assert nx.is_isomorphic(expected_clique, found_clique, edge_match=edge_match), (
            "Clique structure mismatch"
        )
        assert expected_edge == found_edge, "Break edge mismatch"


def test_classify_one_node_no_edges(
    one_node_no_edges_graph: tuple[nx.Graph, str],
) -> None:
    """Check all possible cliques are identified if the initial graph has only one node."""
    graph, break_edge = one_node_no_edges_graph  # Unpack the tuple

    # Finding cliques
    initial_cliques = classify.find_initial_cliques(graph)
    recursive_cliques = classify.find_cliques_recursively(graph)

    # Combine found cliques and expected result
    found_cliques = classify.get_unique_cliques(initial_cliques, recursive_cliques)
    expected_cliques = [(graph, break_edge)]

    # Check the number of identified cliques
    assert len(found_cliques) == 1, "Clique count mismatch"

    # Check that both the graph structures and associated values match
    for (expected_clique, expected_edge), (found_clique, found_edge) in zip(
        expected_cliques, found_cliques, strict=False
    ):
        assert nx.is_isomorphic(expected_clique, found_clique), (
            "Clique structure mismatch"
        )
        assert expected_edge == found_edge, "Clique value mismatch"


def test_classify_two_nodes_no_edges(
    two_nodes_no_edges_graph: nx.Graph, two_nodes_no_edges_cliques: list[nx.Graph]
) -> None:
    """Check all possible cliques are identified if the initial graph has twno nodes and no edges."""
    graph = two_nodes_no_edges_graph

    # Finding cliques
    initial_cliques = (
        classify.find_initial_cliques(graph)
        if len(list(nx.connected_components(graph))) != 1
        else []
    )
    recursive_cliques = classify.find_cliques_recursively(graph)

    # Combine found cliques and expected result
    found_cliques = classify.get_unique_cliques(initial_cliques, recursive_cliques)

    # Check the number of identified cliques
    assert len(found_cliques) == len(two_nodes_no_edges_cliques), (
        "Clique count mismatch"
    )

    # Check that both the graph structures and associated values match
    for (expected_clique, expected_value), (found_clique, found_value) in zip(
        two_nodes_no_edges_cliques, found_cliques, strict=False
    ):
        assert nx.is_isomorphic(expected_clique, found_clique), (
            "Clique structure mismatch"
        )
        assert expected_value == found_value, "Clique value mismatch"


def test_classify_two_nodes_one_edge(
    two_nodes_one_edge_graph: nx.Graph,
    two_nodes_one_edge_cliques: list[tuple[nx.Graph, float | str]],
) -> None:
    """Check all possible cliques are identified if the initial graph has twno nodes and one edges."""
    graph, _break_edge = two_nodes_one_edge_graph

    # Finding cliques
    initial_cliques = classify.find_initial_cliques(graph)
    recursive_cliques = classify.find_cliques_recursively(graph)
    edge_match = iso.numerical_edge_match("coverage", "score")

    found_cliques = classify.get_unique_cliques(initial_cliques, recursive_cliques)

    # Check the number of identified cliques
    assert len(found_cliques) == len(two_nodes_one_edge_cliques), (
        "Clique count mismatch"
    )

    # Check that both the graph structures and associated values match
    for (expected_clique, expected_value), (found_clique, found_value) in zip(
        two_nodes_one_edge_cliques, found_cliques, strict=False
    ):
        assert nx.is_isomorphic(expected_clique, found_clique, edge_match=edge_match), (
            "Clique structure mismatch"
        )
        assert expected_value == found_value, "Clique value mismatch"


def test_classify_complex_graph(
    known_complex_graph: nx.Graph,
    known_complex_cliques: list[tuple[nx.Graph, float]],
) -> None:
    """Check all possible cliques are identified in the initial iteration."""
    graph = known_complex_graph
    initial_cliques = classify.find_initial_cliques(graph)
    recursive_cliques = classify.find_cliques_recursively(graph)
    edge_match = iso.numerical_edge_match("coverage", "score")

    found_cliques = classify.get_unique_cliques(initial_cliques, recursive_cliques)

    # Check the number of identified cliques
    assert len(found_cliques) == len(known_complex_cliques), "Clique count mismatch"

    # Check that both the graph structures and associated values match
    for (expected_clique, expected_value), (found_clique, found_value) in zip(
        known_complex_cliques, found_cliques, strict=False
    ):
        assert nx.is_isomorphic(expected_clique, found_clique, edge_match=edge_match), (
            "Clique structure mismatch"
        )
        assert expected_value == found_value, "Clique value mismatch"
