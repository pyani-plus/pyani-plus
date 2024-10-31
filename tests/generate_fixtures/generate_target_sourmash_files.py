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
"""Generate target files for pyani-plus sourmash tests.

This script can be run with ``./generate_target_sourmash_files.py`` in the script's
directory, or from the project root directory via ``make fixtures``. It will
regenerate and potentially modify test input files under the fixtures directory.

This script generates target files for sourmash comparisons.
Genomes are compared in both directions (forward and reverse)
using `sourmash sketch dna` and `sourmash compare`.

sourmash sketch dna runs with k=31, scaled=1 to find DNA sketches.
Note: By default, sketch dna uses the parameter string "k=31,scaled=1000".
However, these settings fail to estimate ANI values for the current test
set (viral genomes). For testing purposes, we will set the parameters to
"k=31,scaled=1", which does return an estimation of ANI.

sourmash compare runs with --max-containment parameter to compare signatures
and estimate ANI.
"""

# Imports
import subprocess
from itertools import product
from pathlib import Path

from pyani_plus.tools import get_sourmash

# Paths to directories (eg, input sequences, delta and filter)
INPUT_DIR = Path("../fixtures/viral_example")
SIGNATURE_DIR = Path("../fixtures/sourmash/targets/signatures")
COMPARE_DIR = Path("../fixtures/sourmash/targets/compare")

# Running ANIm comparisons (all vs all)
inputs = {_.stem: _ for _ in Path(INPUT_DIR).glob("*.f*")}
comparisons = product(inputs, inputs)

# Cleanup
for file in SIGNATURE_DIR.glob("*.sig"):
    file.unlink()
for file in COMPARE_DIR.glob("*.csv"):
    file.unlink()

sourmash = get_sourmash()
print(f"Using nucmer {sourmash.version} at {sourmash.exe_path}")

for genomes in comparisons:
    stem = "_vs_".join(genomes)
    subprocess.run(
        [
            sourmash.exe_path,
            "sketch",
            "dna",
            "-p",
            "k=31,scaled=1",
            inputs[genomes[1]],
            inputs[genomes[0]],
            "-o",
            SIGNATURE_DIR / (stem + ".sig"),
        ],
        check=True,
    )

    subprocess.run(
        [
            sourmash.exe_path,
            "compare",
            SIGNATURE_DIR / (stem + ".sig"),
            "--csv",
            COMPARE_DIR / (stem + ".csv"),
            "--estimate-ani",
            "--max-containment",
        ],
        check=True,
    )
