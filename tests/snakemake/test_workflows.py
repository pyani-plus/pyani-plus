# The MIT License
#
# Copyright (c) 2024-2025 University of Strathclyde
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
"""Tests for the pyani_plus/snakemake/snakemake_scheduler.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

import re
from pathlib import Path

import pytest

from pyani_plus import setup_logger
from pyani_plus.workflows import (
    ShowProgress,
    ToolExecutor,
    run_snakemake_with_progress_bar,
)


def test_progress_bar_error() -> None:
    """Verify expected error message without database or run-in."""
    msg = "Both database and run_id are required with display as progress bar"
    logger = setup_logger(None)
    with pytest.raises(SystemExit, match=re.escape(msg)):
        run_snakemake_with_progress_bar(
            logger,
            ToolExecutor.slurm,
            "some_workflow.smk",
            [Path("/mnt/shared/answer.tsv")],
            Path(),
            Path(),
            display=ShowProgress.bar,
        )
