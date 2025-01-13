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
"""Generate target files for pyani-plus ANIb (blast) tests.

This script can be run with ``./generate_anib_blast_files.py`` in the script's
directory, or from the project root directory via ``make fixtures``. It will
regenerate and potentially modify test input files under the fixtures directory.

This script generates BLAST databases (from which we take the *.njs files for
test output checking) and then runs blastn using the previously generated
fragmented FASTA files as queries.
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from pyani_plus.methods.anib import BLAST_COLUMNS
from pyani_plus.tools import get_blastn, get_makeblastdb

# Paths to directories (eg, input sequences, delta and filter)
SUBJECT_DIR = Path("../fixtures/viral_example")
QUERY_DIR = NJS_DIR = BLASTN_DIR = Path("../fixtures/viral_example/intermediates/ANIb/")

# Remove pre-existing fixtures before regenerating new ones.
# This is to help with if and when we change the
# example sequences being used.
for file in NJS_DIR.glob("*.njs"):
    file.unlink()

blastn = get_blastn()
makeblastdb = get_makeblastdb()
if blastn.version != makeblastdb.version:
    sys.exit(
        "ERROR - Inconsistent versions, blastn {blastn.version} vs makeblastn {makeblastdb.version}"
    )
print(
    f"Using NCBI BLAST+ {makeblastdb.version} at {blastn.exe_path} and {makeblastdb.exe_path}"
)

count = 0
# Note flexible on input *.fna vs *.fa vs *.fasta, but fixed on output
for subject in SUBJECT_DIR.glob("*.f*"):
    output = NJS_DIR / (subject.stem + ".njs")
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(
            [
                str(makeblastdb.exe_path),
                "-in",
                subject,
                "-title",
                subject.stem,
                "-dbtype",
                "nucl",
                "-out",
                tmp + "/" + subject.stem,
            ],
            check=True,
        )
        shutil.move(
            tmp + "/" + subject.stem + ".njs", NJS_DIR / (subject.stem + ".njs")
        )
        print(f"Collected {subject.stem} BLAST nucleotide database JSON file")
        for query in QUERY_DIR.glob("*-fragments.fna"):
            subprocess.run(
                [
                    str(blastn.exe_path),
                    "-query",
                    query,
                    "-db",
                    tmp + "/" + subject.stem,
                    "-outfmt",
                    # We do NOT use the standard 12 columns
                    "6 " + " ".join(BLAST_COLUMNS),
                    "-out",
                    BLASTN_DIR / f"{query.stem[:-10]}_vs_{subject.stem}.tsv",
                    "-xdrop_gap_final",
                    "150",
                    "-dust",
                    "no",
                    "-evalue",
                    "1e-15",
                    "-task",
                    "blastn",
                ],
                check=True,
            )
        count += 1
print(f"Collected {count} BLAST nucleotide database JSON files")
