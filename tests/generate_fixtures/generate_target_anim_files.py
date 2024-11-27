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
"""Generate target files for pyani-plus ANIm tests.

This script can be run with ``./generate_target_anim_files.py`` in the script's
directory, or from the project root directory via ``make fixtures``. It will
regenerate and potentially modify test input files under the fixtures directory.

This script generates target files for anim comparisons.
Genomes are compared in both directions (forward and reverse)
using nucmer and delta-filter.

nucmer runs with the --mum parameter to find initial maximal unique matches.
delta-filter runs with the -1 parameter to filter only 1-to-1 matches
"""

# Imports
import subprocess
from itertools import product
from pathlib import Path

from pyani_plus.tools import get_delta_filter, get_nucmer

# Generating fixtures for viral example
# Paths to directories (eg, input sequences, delta and filter)
VIRAL_INPUT_DIR = Path("../fixtures/viral_example")
VIRAL_DELTA_DIR = VIRAL_FILTER_DIR = Path(
    "../fixtures/viral_example/intermediates/ANIm"
)

# Running ANIm comparisons (all vs all)
inputs = {_.stem: _ for _ in Path(VIRAL_INPUT_DIR).glob("*.f*")}
comparisons = product(inputs, inputs)

# Cleanup
for file in VIRAL_DELTA_DIR.glob("*.delta"):
    file.unlink()
for file in VIRAL_FILTER_DIR.glob("*.filter"):
    file.unlink()

nucmer = get_nucmer()
delta_filter = get_delta_filter()
print(f"Using nucmer {nucmer.version} at {nucmer.exe_path}")

for genomes in comparisons:
    stem = "_vs_".join(genomes)
    subprocess.run(
        [
            nucmer.exe_path,
            "-p",
            VIRAL_DELTA_DIR / stem,
            "--mum",
            inputs[genomes[1]],
            inputs[genomes[0]],
        ],
        check=True,
    )

    # To redirect using subprocess.run, we need to open the output file and
    # pipe from within the call to stdout
    with (VIRAL_FILTER_DIR / (stem + ".filter")).open("w") as ofh:
        subprocess.run(
            [delta_filter.exe_path, "-1", VIRAL_DELTA_DIR / (stem + ".delta")],
            check=True,
            stdout=ofh,
        )

# Generating fixtures for bad alignments example
# Paths to directories (eg, input sequences, delta and filter)
BAD_ALIGNMENTS_INPUT_DIR = Path("../fixtures/bad_alignments")
BAD_ALIGNMENTS_DELTA_DIR = BAD_ALIGNMENTS_FILTER_DIR = Path(
    "../fixtures/bad_alignments/intermediates/ANIm"
)

# Generating fixtures for genomes that return no alignments
# Running ANIm comparisons.
# Skipping self-to-self comparisons for now.
bad_inputs = {_.stem: _ for _ in Path(BAD_ALIGNMENTS_INPUT_DIR).glob("*.f*")}
bad_comparisons = [
    (genome_1, genome_2)
    for genome_1, genome_2 in product(bad_inputs, bad_inputs)
    if genome_1 != genome_2
]

# Cleanup
for file in BAD_ALIGNMENTS_DELTA_DIR.glob("*.delta"):
    file.unlink()
for file in BAD_ALIGNMENTS_FILTER_DIR.glob("*.filter"):
    file.unlink()

for genomes in bad_comparisons:
    stem = "_vs_".join(genomes)
    subprocess.run(
        [
            nucmer.exe_path,
            "-p",
            BAD_ALIGNMENTS_DELTA_DIR / stem,
            "--mum",
            bad_inputs[genomes[1]],
            bad_inputs[genomes[0]],
        ],
        check=True,
    )

    # To redirect using subprocess.run, we need to open the output file and
    # pipe from within the call to stdout
    with (BAD_ALIGNMENTS_FILTER_DIR / (stem + ".filter")).open("w") as ofh:
        subprocess.run(
            [delta_filter.exe_path, "-1", BAD_ALIGNMENTS_DELTA_DIR / (stem + ".delta")],
            check=True,
            stdout=ofh,
        )
