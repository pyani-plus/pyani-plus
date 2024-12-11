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

This script can be run with
``./generate_target_sourmash_files.py <path_to_inputs_dir> <path_to_output_dir>``
in the script's directory, or from the project root directory via ``make fixtures``.
It will regenerate and potentially modify test input files under thefixtures
directory.

This script generates target files for sourmash comparisons.
Genomes are compared in both directions (forward and reverse)
using `sourmash sketch dna` and `sourmash compare`.

sourmash sketch dna runs with k=31, scaled=300 to find DNA sketches.
Note: By default, sketch dna uses the parameter string "k=31,scaled=1000".
However, these settings fail to estimate ANI values for the current test
set (viral genomes). For testing purposes, we will set the parameters to
"k=31,scaled=300", which does return an estimation of ANI.

sourmash compare runs with --containment parameter to compare signatures
and estimate ANI (we infer the --max-containment value).
"""

# Imports
import subprocess
import sys
from itertools import product
from pathlib import Path

from pyani_plus.tools import get_sourmash

# Paths to directories (eg, input sequences, delta and filter)
INPUT_DIR, OUT_DIR = Path(sys.argv[1]), Path(sys.argv[2])

# Running ANIm comparisons (all vs all)
inputs = {_.stem: _ for _ in INPUT_DIR.glob("*.f*")}
comparisons = product(inputs, inputs)

# Cleanup
for file in OUT_DIR.glob("*.sig"):
    file.unlink()
for file in OUT_DIR.glob("*.csv"):
    file.unlink()

sourmash = get_sourmash()
print(f"Using nucmer {sourmash.version} at {sourmash.exe_path}")

for genome in inputs.values():
    subprocess.run(
        [
            sourmash.exe_path,
            "sketch",
            "dna",
            "-p",
            "k=31,scaled=300",
            genome.resolve().relative_to(OUT_DIR.resolve(), walk_up=True),
            # Seems branchwater requires the name field, which by default is set to
            # the filename (without path) by sourmash scripts singlesketch
            "--name",
            genome.name,
            "-o",
            genome.stem + ".sig",
        ],
        cwd=str(OUT_DIR),  # so paths in .sig file are relative
        check=True,
    )

# This is the traditional sourmash compare output,
# which is used to generate the expected matrices
subprocess.run(
    [
        sourmash.exe_path,
        "compare",
        *sorted(OUT_DIR.glob("*.sig")),
        "--csv",
        OUT_DIR / "sourmash.csv",
        "--estimate-ani",
        "--containment",
        "-k=31",
    ],
    check=True,
)
