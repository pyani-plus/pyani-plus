#!/usr/bin/env python3
#
# The MIT License
#
# Copyright (c) 2024-2026 University of Strathclyde
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
"""Generate target files for pyani-plus skani tests.

This script can be run with
``python generate_target_skani_files.py <path_to_inputs_dir> <path_to_output_dir>``
in the script's directory, or from the project root directory
via ``make fixtures``. It will regenerate and potentially modify test
input files under the fixtures directory.

This script generates target files for skani comparisons.
Genomes are compared in both directions (forward and reverse)
using skani dist
"""

# Imports
import subprocess
import sys
from itertools import product
from pathlib import Path

from pyani_plus.tools import get_skani
from pyani_plus.utils import file_md5sum

# Generating fixtures
# Assign variable to specify Paths to directories (eg. input, output etc.)
INPUT_DIR, OUT_DIR = Path(sys.argv[1]), Path(sys.argv[2])
OUT_DIR.mkdir(exist_ok=True, parents=True)

# Running ANIm comparisons (all vs all)
inputs = {_.stem: _ for _ in INPUT_DIR.glob("*.f*")}
comparisons = product(inputs, inputs)

skani = get_skani()
print(f"Using skani {skani.version} at {skani.exe_path}")

# skani writes to stdout.
# To redirect using subprocess.run, we need to open the output file and
# pipe from within the call to stdout
for genomes in comparisons:
    stem = "_vs_".join(file_md5sum(inputs[_]) for _ in genomes)
    with (OUT_DIR / (stem + ".skani")).open("w") as ofh:
        subprocess.run(
            [
                skani.exe_path,
                "dist",
                inputs[genomes[1]],
                inputs[genomes[0]],
            ],
            check=True,
            stdout=ofh,
        )
