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
"""Test methods for calculating ANIm

These tests are intended to be run from the repository root using:

make test
"""

# Required to support pytest automated testing
import pytest

import pandas as pd
from pathlib import Path

from pyani_plus.methods import method_anim


@pytest.fixture
def anim_results(dir_anim_results):
    """Expected results for anim."""
    return pd.read_csv(dir_anim_results / "deltadir_result.csv", index_col=0)


def compare_anim_results(filterfile, datadir):
    """Compare delta-filter files and expected anim results.

    This function expects delta-filter file and dataframe with expected
    anim results as input and returns True if the anim values is the same,
    and False if the two files differ.
    """

    reference, query = filterfile.stem.split("_vs_")

    if round(method_anim.parse_delta(filterfile)[2], 6) != datadir.at[reference, query]:
        return False

    return True


def test_anim_parsing(anim_nucmer_targets_filter, anim_results):
    """Check parsing of test NUCmer .filter file."""

    for fname in anim_nucmer_targets_filter:
        assert compare_anim_results(fname, anim_results)
