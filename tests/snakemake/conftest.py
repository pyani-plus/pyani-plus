# Copyright (c) 2024-present University of Strathclyde
# Author: Peter Cock
#
# Contact:
# peter.cock@strath.ac.uk
#
# Peter Cock,
# Strathclyde Institute of Pharmaceutical and Biomedical Sciences
# The University of Strathclyde
# 161 Cathedral Street
# Glasgow
# G4 0RE
# Scotland,
# UK
#
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
"""Pytest configuration file for our snakemake tests."""

import os

import pytest


@pytest.fixture
def snakemake_cores() -> int:
    """How many CPU cores/threads for snakemake, 8 unless capped by those available to use.

    Returns 8 cores, unless capped by the system limits (e.g. SLURM might only give 4 cores,
    or GitHub Actions only 2 cores).
    """
    try:
        # This will take into account SLURM limits,
        # so don't need to check $SLURM_CPUS_PER_TASK explicitly.
        # Probably don't need to check $NSLOTS on SGE either.
        available = len(os.sched_getaffinity(0))  # type: ignore[attr-defined]
    except AttributeError:
        # Unavailable on macOS or Windows, use this instead
        # Can return None
        cpus = os.cpu_count()
        if not cpus:
            msg = "Cannot determine CPU count"
            raise RuntimeError(msg) from None
        available = cpus
    return min(8, available)
