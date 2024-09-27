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
"""Generate target files for pyani-plus ANIb (blast) tests.

This script can be run with ``./generate_anib_blastdb_njs_files.py`` in the script's
directory, or from the project root directory via ``make fixtures``. It will
regenerate and potentially modify test input files under the fixtures directory.

This script generates fragmented FASTA files which are used as the query files
with NCBI BLAST+ command blastn against database of the unfragmented FASTA files.
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

from pyani_plus.tools import get_makeblastdb

# Paths to directories (eg, input sequences, delta and filter)
INPUT_DIR = Path("../fixtures/viral_example")
NJS_DIR = Path("../fixtures/anib/blastdb")

makeblastdb = get_makeblastdb()
print(f"Using NCBI BLAST+ {makeblastdb.version} at {makeblastdb.exe_path}")  # noqa: T201

count = 0
# Note flexible on input *.fna vs *.fa vs *.fasta, but fixed on output
for fasta in INPUT_DIR.glob("*.f*"):
    output = NJS_DIR / (fasta.stem + ".njs")
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(
            [
                makeblastdb.exe_path,
                "-in",
                fasta,
                "-title",
                fasta.stem,
                "-dbtype",
                "nucl",
                "-out",
                tmp + "/" + fasta.stem,
            ],
            check=True,
        )
        shutil.move(tmp + "/" + fasta.stem + ".njs", NJS_DIR / (fasta.stem + ".njs"))
        count += 1
    print(f"Collected {fasta.stem} BLAST nucleotide database JSON file")  # noqa: T201
print(f"Collected {count} BLAST nucleotide database JSON files")  # noqa: T201
