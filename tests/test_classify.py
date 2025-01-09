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
import pytest

from pyani_plus import classify


@pytest.fixture
def expected_complete_graph() -> nx.Graph:
    """Example of complete graph (viral example)."""  # noqa: D401
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
            "identity": 0.9946424059000001,
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
def expected_cluster_info() -> classify.ClusterInfo:
    """Example of expected classify ClusterInfo."""  # noqa: D401
    return classify.ClusterInfo(
        n_nodes=3,
        members=sorted(["OP073605", "MGV-GENOME-0264574", "MGV-GENOME-0266457"]),
        min_cov=0.6774176803,
        min_identity=0.9946424059000001,
        clique=True,
        singleton=False,
    )


def test_construct_complete_graph(
    input_genomes_tiny: Path, expected_complete_graph: nx.Graph
) -> None:
    """Check construction of complete graph."""
    # Comparison function for a numerical edge attribute.
    edge_match = iso.numerical_edge_match("coverage", "identity")

    # Check the isomorphism of a graph with the edge_match function
    assert nx.is_isomorphic(
        expected_complete_graph,
        classify.construct_complete_graph(input_genomes_tiny / "viral_database.db", 1),
        edge_match=edge_match,
    )


def test_edges_below_threshold(
    input_genomes_tiny: Path, expected_complete_graph: nx.Graph
) -> None:
    """Check removal of edges below specified threshold."""
    edge_match = iso.numerical_edge_match("coverage", "identity")
    complete_graph = classify.construct_complete_graph(
        input_genomes_tiny / "viral_database.db", 1
    )

    # No edges have a weight below 50% coverage, so none should be removed,
    # ensuring the final graph matches the expected complete graph.
    assert nx.is_isomorphic(
        expected_complete_graph,
        classify.remove_edges_below_threshold(complete_graph, "coverage", 0.5),
        edge_match=edge_match,
    )


def test_analyse_subgraphs(
    input_genomes_tiny: Path, expected_cluster_info: classify.ClusterInfo
) -> None:
    """Check classification of subgraphs in a graph."""
    complete_graph = classify.construct_complete_graph(
        input_genomes_tiny / "viral_database.db", 1
    )
    assert expected_cluster_info == classify.analyse_subgraphs(complete_graph)[0]
