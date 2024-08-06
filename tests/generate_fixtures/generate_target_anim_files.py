# (c) University of Strathclyde 2024-present
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
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
"""Generate target files for pyani-plus ANIm tests.

This bash script generates target files for anim comparisons.
Genomes are compared in both directions (forward and reverse)
using nucmer and delta-filter.

nucmer runs with the --mum parameter to find initial maximal unique matches.
delta-filter runs with the -1 parameter to filter only 1-to-1 matches
"""

# Imports
import subprocess
from itertools import permutations
from pathlib import Path

# Paths to directories (eg, input sequences, delta and filter)
INPUT_DIR = Path("../fixtures/sequences")
DELTA_DIR = Path("../fixtures/anim/targets/delta")
FILTER_DIR = Path("../fixtures/anim/targets/filter")

# Running ANIm comparisons
comparisons = permutations([_.stem for _ in Path(INPUT_DIR).glob("*")], 2)
inputs = {_.stem: _ for _ in Path(INPUT_DIR).glob("*")}

for genomes in comparisons:
    stem = "_vs_".join(genomes)
    subprocess.run(
        [  # noqa: S607
            "nucmer",
            "-p",
            DELTA_DIR / stem,
            "--mum",
            inputs[genomes[0]],
            inputs[genomes[1]],
        ],
        check=True,
    )

    # To redirect using subprocess.run, we need to open the output file and
    # pipe from within the call to stdout
    with (FILTER_DIR / (stem + ".filter")).open("w") as ofh:
        subprocess.run(
            ["delta-filter", "-1", DELTA_DIR / (stem + ".delta")],  # noqa: S607
            check=True,
            stdout=ofh,
        )
