#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# The MIT License
#
# Copyright (c) University of Strathclyde
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

"""Test run_anim_workflow.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

#Set Up 
from pathlib import Path
import sys
sys.path.append("..")
from itertools import permutations
from scripts.run_anim_workflow import get_target_files
import pytest


config_args = {
        "outdir": "../issue_1/output",
        "indir": "../issue_1/input",
        "cores": 8,
        "mode": 'maxmatch'
    }

 


@pytest.fixture
def target_files():
    path_to_target_files = Path("test_targets/anim")
    return [fname.name for fname in sorted(list(path_to_target_files.glob("*"))) if fname.suffix == '.filter']




# Test function that uses the fixture
def test_get_target_files(target_files):
    path_to_inputs = Path("test_input/anim")

    assert target_files ==  get_target_files(path_to_inputs)
    # Add more assertions or tests as needed



