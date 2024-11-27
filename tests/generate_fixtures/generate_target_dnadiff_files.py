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
"""Generate target files for pyani-plus dnadiff tests.

This script can be run with ``./generate_target_dnadiff_files.py`` in the script's
directory, or from the project root directory via ``make fixtures``. It will
regenerate and potentially modify test input files under the fixtures directory.

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
import shutil
import subprocess
import tempfile
from itertools import product
from pathlib import Path

from pyani_plus.tools import (
    get_delta_filter,
    get_dnadiff,
    get_nucmer,
    get_show_coords,
    get_show_diff,
)

# Generating fixtures for viral example
# Paths to directories (eg. input sequences, outputs for delta, filter...)
INPUT_DIR = Path("../fixtures/viral_example")
DELTA_DIR = FILTER_DIR = SHOW_DIFF_DIR = SHOW_COORDS_DIR = DNADIFF_DIR = Path(
    "../fixtures/viral_example/intermediates/dnadiff/"
)

# Running comparisons (all vs all)
inputs = {_.stem: _ for _ in sorted(Path(INPUT_DIR).glob("*.f*"))}
comparisons = product(inputs, inputs)

# # Cleanup
# This is to help with if and when we change the
# example sequences being used.
for file in DELTA_DIR.glob("*.delta"):
    file.unlink()
for file in FILTER_DIR.glob("*.filter"):
    file.unlink()
for file in SHOW_DIFF_DIR.glob("*.qdiff"):
    file.unlink()
for file in SHOW_COORDS_DIR.glob("*.mcoords"):
    file.unlink()
for file in DNADIFF_DIR.glob("*.report"):
    file.unlink()

nucmer = get_nucmer()
delta_filter = get_delta_filter()
show_coords = get_show_coords()
show_diff = get_show_diff()
dnadiff = get_dnadiff()
print(f"Using nucmer {nucmer.version} at {nucmer.exe_path}")
print(f"Using dnadiff {dnadiff.version} at {dnadiff.exe_path}")

for genomes in comparisons:
    stem = "_vs_".join(genomes)
    subprocess.run(
        [
            nucmer.exe_path,
            "-p",
            DELTA_DIR / stem,
            "--maxmatch",
            inputs[genomes[1]],
            inputs[genomes[0]],
        ],
        check=True,
    )

    # To redirect using subprocess.run, we need to open the output file and
    # pipe from within the call to stdout
    with (FILTER_DIR / (stem + ".filter")).open("w") as ofh:
        subprocess.run(
            [delta_filter.exe_path, "-m", DELTA_DIR / (stem + ".delta")],
            check=True,
            stdout=ofh,
        )

    with (SHOW_DIFF_DIR / (stem + ".qdiff")).open("w") as ofh:
        subprocess.run(
            [show_diff.exe_path, "-qH", FILTER_DIR / (stem + ".filter")],
            check=True,
            stdout=ofh,
        )

    with (SHOW_COORDS_DIR / (stem + ".mcoords")).open("w") as ofh:
        subprocess.run(
            [show_coords.exe_path, "-rclTH", FILTER_DIR / (stem + ".filter")],
            check=True,
            stdout=ofh,
        )
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(
            [
                dnadiff.exe_path,
                "-p",
                tmp + "/" + stem,
                inputs[genomes[1]],
                inputs[genomes[0]],
            ],
            check=True,
        )
        shutil.move(tmp + "/" + stem + ".report", DNADIFF_DIR / (stem + ".report"))


# Generating fixtures for bad alignments example
# Paths to directories (eg. input sequences, outputs for delta, filter...)
BAD_ALIGNMENTS_INPUT_DIR = Path("../fixtures/bad_alignments")
BAD_ALIGNMENTS_DELTA_DIR = BAD_ALIGNMENTS_FILTER_DIR = BAD_ALIGNMENTS_SHOW_DIFF_DIR = (
    BAD_ALIGNMENTS_SHOW_COORDS_DIR
) = BAD_ALIGNMENTS_DNADIFF_DIR = Path(
    "../fixtures/bad_alignments/intermediates/dnadiff/"
)

# Generating fixtures for genomes that return no alignments
# Running dnadiff comparisons.
# Skipping self-to-self comparisons for now.
bad_inputs = {_.stem: _ for _ in sorted(Path(BAD_ALIGNMENTS_INPUT_DIR).glob("*.f*"))}
bad_comparisons = [
    (genome_1, genome_2)
    for genome_1, genome_2 in product(bad_inputs, bad_inputs)
    if genome_1 != genome_2
]

# Cleanup
# This is to help with if and when we change the
# example sequences being used.
for file in BAD_ALIGNMENTS_DELTA_DIR.glob("*.delta"):
    file.unlink()
for file in BAD_ALIGNMENTS_FILTER_DIR.glob("*.filter"):
    file.unlink()
for file in BAD_ALIGNMENTS_SHOW_DIFF_DIR.glob("*.qdiff"):
    file.unlink()
for file in BAD_ALIGNMENTS_SHOW_COORDS_DIR.glob("*.mcoords"):
    file.unlink()
for file in BAD_ALIGNMENTS_DNADIFF_DIR.glob("*.report"):
    file.unlink()

nucmer = get_nucmer()
delta_filter = get_delta_filter()
show_coords = get_show_coords()
show_diff = get_show_diff()
dnadiff = get_dnadiff()
print(f"Using nucmer {nucmer.version} at {nucmer.exe_path}")
print(f"Using dnadiff {dnadiff.version} at {dnadiff.exe_path}")

for genomes in bad_comparisons:
    stem = "_vs_".join(genomes)
    subprocess.run(
        [
            nucmer.exe_path,
            "-p",
            BAD_ALIGNMENTS_DELTA_DIR / stem,
            "--maxmatch",
            bad_inputs[genomes[1]],
            bad_inputs[genomes[0]],
        ],
        check=True,
    )

    # To redirect using subprocess.run, we need to open the output file and
    # pipe from within the call to stdout
    with (BAD_ALIGNMENTS_FILTER_DIR / (stem + ".filter")).open("w") as ofh:
        subprocess.run(
            [delta_filter.exe_path, "-m", BAD_ALIGNMENTS_DELTA_DIR / (stem + ".delta")],
            check=True,
            stdout=ofh,
        )

    with (BAD_ALIGNMENTS_SHOW_DIFF_DIR / (stem + ".qdiff")).open("w") as ofh:
        subprocess.run(
            [show_diff.exe_path, "-qH", BAD_ALIGNMENTS_FILTER_DIR / (stem + ".filter")],
            check=True,
            stdout=ofh,
        )

    with (BAD_ALIGNMENTS_SHOW_COORDS_DIR / (stem + ".mcoords")).open("w") as ofh:
        subprocess.run(
            [
                show_coords.exe_path,
                "-rclTH",
                BAD_ALIGNMENTS_FILTER_DIR / (stem + ".filter"),
            ],
            check=True,
            stdout=ofh,
        )
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(
            [
                dnadiff.exe_path,
                "-p",
                tmp + "/" + stem,
                bad_inputs[genomes[1]],
                bad_inputs[genomes[0]],
            ],
            check=True,
        )
        shutil.move(
            tmp + "/" + stem + ".report",
            BAD_ALIGNMENTS_DNADIFF_DIR / (stem + ".report"),
        )
