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
"""Generate target files for pyani-plus lz-ani tests.

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
import tempfile
from itertools import product
from pathlib import Path

from pyani_plus.tools import get_lzani
from pyani_plus.utils import file_md5sum

# Generating fixtures
# Assign variable to specify Paths to directories (eg. input, output etc.)
INPUT_DIR, OUT_DIR = Path(sys.argv[1]), Path(sys.argv[2])
OUT_DIR.mkdir(exist_ok=True, parents=True)

# Running ANIm comparisons (all vs all)
inputs = {_.stem: _ for _ in INPUT_DIR.glob("*.f*")}
comparisons = product(inputs, inputs)

lzani = get_lzani()
print(f"Using lzani {lzani.version} at {lzani.exe_path}")

# lz-ani is more complex to organise than most other tools as it does not
# directly allow for pairwise comparisons, and it will by default treat
# individual contigs as separate genomes. To get around this, we provide
# the --in-txt argument with a text file specifying paths to the two
# genomes we want to compare.o stdout
for genomes in comparisons:
    stem = "_vs_".join(file_md5sum(inputs[_]) for _ in genomes)
    # Write a text file with the paths to the two genomes to compare, one per line
    # to a temporary file, which we will then pass to lz-ani with the --in-txt argument
    with tempfile.TemporaryDirectory() as tmpdir:
        in_txt_path = Path(tmpdir) / (stem + "_in.txt")

        with in_txt_path.open("w") as ofh:
            ofh.write(str(inputs[genomes[0]]) + "\n")
            ofh.write(str(inputs[genomes[1]]) + "\n")

        # Run lz-ani comparison
        subprocess.run(
            [
                str(lzani.exe_path),
                "all2all",
                "--multisample-fasta",
                "false",
                "--in-dir",
                str(inputs[genomes[0]].parent),
                "--in-txt",
                str(in_txt_path),
                "-o",
                str(OUT_DIR / (stem + ".tsv")),
            ],
            check=True,
        )
