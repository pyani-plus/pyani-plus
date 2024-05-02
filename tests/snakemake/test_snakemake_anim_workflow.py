#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# The MIT License
#
# Copyright (c) 2024-present University of Strathclyde
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
"""Test snakemake workflow for ANIm

These tests are intended to be run from the repository root using:

make test
"""

import shutil  # We need this for filesystem operations

# Required to support pytest automated testing
import pytest

# We're testing the workflow, as called through the pyani_plus
# wrapper code, so import only that module
from pyani_plus.snakemake import anim


@pytest.fixture
def config_filter_args(anim_nucmer_targets_filter_outdir, input_genomes_small):
    """Configuration settings for testing snakemake filter rule.

    We take the output directories for the MUMmer filter output and the
    small set of input genomes as arguments.
    """
    return {
        "outdir": anim_nucmer_targets_filter_outdir,
        "indir": input_genomes_small,
        "cores": 8,
        "mode": "mum",
    }


@pytest.fixture
def config_delta_args(anim_nucmer_targets_delta_outdir, input_genomes_small):
    """Configuration settings for testing snakemake delta rule.

    We take the output directories for the MUMmer delta output and the
    small set of input genomes as arguments.
    """
    return {
        "outdir": anim_nucmer_targets_delta_outdir,
        "indir": input_genomes_small,
        "cores": 8,
        "mode": "mum",
    }


def compare_files_with_skip(file1, file2, skip=1):
    """Compare two files, line by line, except for the first line.

    This function expects two text files as input and returns True if the content
    of the files is the same, and False if the two files differ.
    """

    with file1.open() as if1, file2.open() as if2:
        for line1, line2 in zip(
            if1.readlines()[skip:], if2.readlines()[skip:], strict=False
        ):
            if line1 != line2:
                return False

    return True


def test_snakemake_rule_filter(
    anim_nucmer_targets_filter,
    anim_nucmer_targets_filter_indir,
    anim_nucmer_targets_filter_outdir,
    config_filter_args,
):
    """Test nucmer filter snakemake wrapper

    Checks that the filter rule in the ANIm snakemake wrapper gives the
    expected output.

    If the output directory exists (i.e. the make clean_tests rule has not
    been run), the tests will automatically pass as snakemake will not
    attempt to re-run the rule. That would prevent us from seeing any
    introduced bugs, so we force re-running the rule by deleting the
    output directory before running the tests.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(anim_nucmer_targets_filter_outdir, ignore_errors=True)

    # Run snakemake wrapper
    anim.run_workflow(anim_nucmer_targets_filter, config_filter_args)

    # Check output against target fixtures
    for fname in anim_nucmer_targets_filter:
        assert compare_files_with_skip(
            anim_nucmer_targets_filter_indir / fname,
            anim_nucmer_targets_filter_outdir / fname,
        )


def test_snakemake_rule_delta(
    anim_nucmer_targets_delta,
    anim_nucmer_targets_delta_indir,
    anim_nucmer_targets_delta_outdir,
    config_delta_args,
):
    """Test nucmer delta snakemake wrapper

    Checks that the delta rule in the ANIm snakemake wrapper gives the
    expected output.

    If the output directory exists (i.e. the make clean_tests rule has not
    been run), the tests will automatically pass as snakemake will not
    attempt to re-run the rule. That would prevent us from seeing any
    introduced bugs, so we force re-running the rule by deleting the
    output directory before running the tests.
    """
    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(anim_nucmer_targets_delta_outdir, ignore_errors=True)

    # Run snakemake wrapper
    anim.run_workflow(anim_nucmer_targets_delta, config_delta_args)

    # Check output against target fixtures
    for fname in anim_nucmer_targets_delta:
        assert compare_files_with_skip(
            anim_nucmer_targets_delta_indir / fname,
            anim_nucmer_targets_delta_outdir / fname,
        )
