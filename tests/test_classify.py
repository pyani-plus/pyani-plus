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

from pathlib import Path

import networkx as nx  # type: ignore  # noqa: PGH003
import networkx.algorithms.isomorphism as iso
import pandas as pd
import pytest

from pyani_plus import classify as method_classify


@pytest.fixture
def expected_complete_graph() -> nx.Graph:
    """Return complete graph (viral example) with average identity and lowest genome coverage."""
    graph = nx.Graph()
    nodes = ["OP073605", "MGV-GENOME-0264574", "MGV-GENOME-0266457"]
    graph.add_nodes_from(nodes)

    edges = [
        ("OP073605", "MGV-GENOME-0264574"),
        ("OP073605", "MGV-GENOME-0266457"),
        ("MGV-GENOME-0264574", "MGV-GENOME-0266457"),
    ]

    edge_attributes = {
        ("OP073605", "MGV-GENOME-0264574"): {
            "identity": 0.999310,
            "coverage": 0.6774176803,
        },
        ("OP073605", "MGV-GENOME-0266457"): {
            "identity": 0.9946424059,
            "coverage": 0.68465,
        },
        ("MGV-GENOME-0264574", "MGV-GENOME-0266457"): {
            "identity": 0.996249,
            "coverage": 0.989443,
        },
    }

    graph.add_edges_from(edges)
    nx.set_edge_attributes(graph, edge_attributes)

    return graph


@pytest.fixture
def expected_cliques() -> nx.Graph:
    """Return list of expected cliques (viral example)."""
    c1 = nx.Graph()
    nodes = ["OP073605", "MGV-GENOME-0264574", "MGV-GENOME-0266457"]
    c1.add_nodes_from(nodes)

    edges = [
        ("OP073605", "MGV-GENOME-0264574"),
        ("OP073605", "MGV-GENOME-0266457"),
        ("MGV-GENOME-0264574", "MGV-GENOME-0266457"),
    ]

    edge_attributes_c1 = {
        ("OP073605", "MGV-GENOME-0264574"): {
            "identity": 0.999310,
            "coverage": 0.6774176803,
        },
        ("OP073605", "MGV-GENOME-0266457"): {
            "identity": 0.9946424059,
            "coverage": 0.68465,
        },
        ("MGV-GENOME-0264574", "MGV-GENOME-0266457"): {
            "identity": 0.996249,
            "coverage": 0.989443,
        },
    }

    c1.add_edges_from(edges)
    nx.set_edge_attributes(c1, edge_attributes_c1)

    c2 = nx.Graph()
    nodes = ["MGV-GENOME-0264574", "MGV-GENOME-0266457"]
    c2.add_nodes_from(nodes)

    edges = [
        ("MGV-GENOME-0264574", "MGV-GENOME-0266457"),
    ]

    edge_attributes_c2 = {
        ("MGV-GENOME-0264574", "MGV-GENOME-0266457"): {
            "identity": 0.996249,
            "coverage": 0.989443,
        },
    }

    c2.add_edges_from(edges)
    nx.set_edge_attributes(c2, edge_attributes_c2)

    return [c1, c2]


def test_construct_complete_graph(
    input_genomes_tiny: Path, expected_complete_graph: nx.Graph
) -> None:
    """Check construction of complete graph."""
    # Matrices
    cov_matrix = pd.read_csv(
        input_genomes_tiny / "matrices/ANIm_coverage.tsv", sep="\t", index_col=0
    )
    id_matrix = pd.read_csv(
        input_genomes_tiny / "matrices/ANIm_identity.tsv", sep="\t", index_col=0
    )
    # Comparison function for a numerical edge attribute.
    edge_match = iso.numerical_edge_match("coverage", "identity")

    # Check the isomorphism of a graph with the edge_match function
    assert nx.is_isomorphic(
        expected_complete_graph,
        method_classify.construct_complete_graph(cov_matrix, id_matrix),
        edge_match=edge_match,
    )


def test_is_clique(expected_complete_graph: nx.Graph) -> None:
    """Check cliques are identified."""
    assert method_classify.is_clique(expected_complete_graph) is True


def test_remove_lowest_edge(expected_complete_graph: nx.Graph) -> None:
    """Check removal of edges threshold."""
    expected_graph_no_edge = expected_complete_graph.copy()
    expected_graph_no_edge.remove_edge("OP073605", "MGV-GENOME-0264574")
    edge_match = iso.numerical_edge_match("coverage", "identity")

    assert nx.is_isomorphic(
        expected_graph_no_edge,
        method_classify.remove_lowest_edge(expected_complete_graph, "coverage"),
        edge_match=edge_match,
    )


def test_find_cliques(
    expected_complete_graph: nx.Graph, expected_cliques: list
) -> None:
    """Check all possible cliques are identified."""
    found_cliques = method_classify.find_cliques(expected_complete_graph, "coverage")
    edge_match = iso.numerical_edge_match("coverage", "identity")

    # Check number of identified cliques
    assert len(expected_cliques) == len(found_cliques), "Clique count mismatch"
    # Check the clique structures
    for expected_clique, found_clique in zip(
        expected_cliques, found_cliques, strict=False
    ):
        assert nx.is_isomorphic(expected_clique, found_clique, edge_match=edge_match), (
            "Clique structure mismatch"
        )
