#!/usr/bin/env python3
#
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
"""Generate target files for pyani-plus dnadiff tests.

This script can be run with
``./generate_target_dnadiff_files.py <path_to_input_dir> <path_to_output_dir>``
in the script's directory, or from the project root directory via ``make fixtures``.
It will regenerate and potentially modify test input files under the fixtures directory.

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
import sys
from itertools import product
from pathlib import Path

from pyani_plus.tools import (
    get_delta_filter,
    get_nucmer,
    get_show_coords,
    get_show_diff,
)
from pyani_plus.utils import file_md5sum

# Generating fixtures for viral example
# Paths to directories (eg. input sequences, outputs for delta, filter...)
INPUT_DIR, OUT_DIR = Path(sys.argv[1]), Path(sys.argv[2])

if not OUT_DIR.is_dir():
    OUT_DIR.mkdir()

# Running comparisons (all vs all)
inputs = {_.stem: _ for _ in sorted(Path(INPUT_DIR).glob("*.f*"))}
stem_hashes = {file.stem: file_md5sum(file) for file in INPUT_DIR.glob("*.f*")}
comparisons = product(inputs, inputs)

# # Cleanup
# This is to help with if and when we change the
# example sequences being used.
for file in OUT_DIR.glob("*.delta"):
    file.unlink()
for file in OUT_DIR.glob("*.filter"):
    file.unlink()
for file in OUT_DIR.glob("*.qdiff"):
    file.unlink()
for file in OUT_DIR.glob("*.mcoords"):
    file.unlink()
for file in OUT_DIR.glob("*.report"):
    file.unlink()

nucmer = get_nucmer()
delta_filter = get_delta_filter()
show_coords = get_show_coords()
show_diff = get_show_diff()
print(f"Using nucmer {nucmer.version} at {nucmer.exe_path}")

for genomes in comparisons:
    stem = "_vs_".join(stem_hashes[_] for _ in genomes)
    subprocess.run(
        [
            nucmer.exe_path,
            "-p",
            OUT_DIR / stem,
            "--maxmatch",
            inputs[genomes[1]],
            inputs[genomes[0]],
        ],
        check=True,
    )

    # To redirect using subprocess.run, we need to open the output file and
    # pipe from within the call to stdout
    with (OUT_DIR / (stem + ".filter")).open("w") as ofh:
        subprocess.run(
            [delta_filter.exe_path, "-m", OUT_DIR / (stem + ".delta")],
            check=True,
            stdout=ofh,
        )

    with (OUT_DIR / (stem + ".qdiff")).open("w") as ofh:
        subprocess.run(
            [show_diff.exe_path, "-qH", OUT_DIR / (stem + ".filter")],
            check=True,
            stdout=ofh,
        )

    with (OUT_DIR / (stem + ".mcoords")).open("w") as ofh:
        subprocess.run(
            [show_coords.exe_path, "-rclTH", OUT_DIR / (stem + ".filter")],
            check=True,
            stdout=ofh,
        )
