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
"""Generate target files for pyani-plus dnadiff tests.

nucmer runs with the --maxmatch parameter to find initial matches
regardless of their uniqueness.
delta-filter runs with the -m parameter to filter only M-to-M
matches.
show-diff runs with the -q parameter to obtain structural
differences for each sequence in the reference.

NOTE: Depending on what values we wish to report using dnadiff
subcommand, we might need to generate show-diff files for all
sequences in the query too.
"""

# Imports
import subprocess
from itertools import permutations
from pathlib import Path

# Paths to directories (eg. input sequences, outputs for delta, filter...)
INPUT_DIR = Path("../fixtures/sequences")
DELTA_DIR = Path("../fixtures/dnadiff/targets/delta")
FILTER_DIR = Path("../fixtures/dnadiff/targets/filter")
SHOW_DIFF_DIR = Path("../fixtures/dnadiff/targets/show_diff")
SHOW_COORDS_DIR = Path("../fixtures/dnadiff/targets/show_coords")

# Running comparisons
comparisons = permutations([_.stem for _ in Path(INPUT_DIR).glob("*")], 2)
inputs = {_.stem: _ for _ in Path(INPUT_DIR).glob("*")}

for genomes in comparisons:
    stem = "_vs_".join(genomes)
    subprocess.run(
        [  # noqa: S607
            "nucmer",
            "-p",
            DELTA_DIR / stem,
            "--maxmatch",
            inputs[genomes[0]],
            inputs[genomes[1]],
        ],
        check=True,
    )

    # To redirect using subprocess.run, we need to open the output file and
    # pipe from within the call to stdout
    with (FILTER_DIR / (stem + ".filter")).open("w") as ofh:
        subprocess.run(
            ["delta-filter", "-m", DELTA_DIR / (stem + ".delta")],  # noqa: S607
            check=True,
            stdout=ofh,
        )

    with (SHOW_DIFF_DIR / (stem + ".rdiff")).open("w") as ofh:
        subprocess.run(
            ["show-diff", "-rH", FILTER_DIR / (stem + ".filter")],  # noqa: S607
            check=True,
            stdout=ofh,
        )

    with (SHOW_COORDS_DIR / (stem + ".mcoords")).open("w") as ofh:
        subprocess.run(
            ["show-coords", FILTER_DIR / (stem + ".filter")],  # noqa: S607
            check=True,
            stdout=ofh,
        )
