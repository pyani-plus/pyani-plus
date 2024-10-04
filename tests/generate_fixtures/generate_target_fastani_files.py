#!/usr/bin/env python3
#
# The MIT License
#
# Copyright (c) 2024 University of Strathclyde
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
"""Generate target files for pyani-plus fastANI tests.

This script can be run with ``./generate_target_fastani_files.py`` in the script's
directory, or from the project root directory via ``make fixtures``. It will
regenerate and potentially modify test input files under the fixtures directory.

This script generates target files for fastani comparisons.
Genomes are compared in both directions (forward and reverse)
using fastANI.
"""

# Imports
import subprocess
from itertools import product
from pathlib import Path

from pyani_plus.tools import get_fastani

# Parameters (eg, input sequences, fastANI outputs, k-mer sizes...)
INPUT_DIR = Path("../fixtures/viral_example")
FASTANI_DIR = Path("../fixtures/fastani/targets")
FRAG_LEN = 3000
KMER_SIZE = 16
MIN_FRAC = 0.2

# Running comparisons
inputs = {_.stem: _ for _ in Path(INPUT_DIR).glob("*.f*")}
comparisons = product(inputs, inputs)

fastani = get_fastani()
print(f"Using fastANI {fastani.version} at {fastani.exe_path}")  # noqa: T201

for genomes in comparisons:
    stem = "_vs_".join(genomes)
    subprocess.run(
        [
            fastani.exe_path,
            "-q",
            inputs[genomes[0]],
            "-r",
            inputs[genomes[1]],
            "-o",
            FASTANI_DIR / (stem + ".fastani"),
            "--fragLen",
            str(FRAG_LEN),
            "-k",
            str(KMER_SIZE),
            "--minFraction",
            str(MIN_FRAC),
        ],
        check=True,
    )
