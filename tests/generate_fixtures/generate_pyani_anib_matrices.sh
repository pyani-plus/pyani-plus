#!/bin/bash
set -euo pipefail
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

echo "This is intended to be used with pyANI version 0.2.13.1, we have:"
average_nucleotide_identity.py --version
# Abort if not matched (via set -e above)
average_nucleotide_identity.py --version | grep "pyani 0\.2\." > /dev/null

# Make a temp subdir
rm -rf tmp_pyani_anib
mkdir tmp_pyani_anib
cd tmp_pyani_anib

echo "Setting up input genomes"
for FASTA in ../../../tests/fixtures/viral_example/*.f*; do
    ln -s "$FASTA" "${FASTA##*/}"
done

echo "Run ANIb comparisions..."
mkdir output
average_nucleotide_identity.py -m ANIb -i . -o output -v -l ANIb.log --force # --labels labels.txt --classes classes.txt

echo "Collecting output for test fixtures..."
# The output names here are based on pyANI v0.3 conventions, the pyani-plus names are shorter
# Note the legacy pyANI output is not sorted by filename stem
mv output/ANIb_percentage_identity.tab ../../fixtures/viral_example/matrices/ANIb_identity.tsv
mv output/ANIb_alignment_coverage.tab ../../fixtures/viral_example/matrices/ANIb_coverage.tsv
mv output/ANIb_alignment_lengths.tab ../../fixtures/viral_example/matrices/ANIb_aln_lengths.tsv
mv output/ANIb_hadamard.tab ../../fixtures/viral_example/matrices/ANIb_hadamard.tsv
mv output/ANIb_similarity_errors.tab ../../fixtures/viral_example/matrices/ANIb_sim_errors.tsv

#Remove temp subdir
cd ..
rm -rf tmp_pyani_anib
echo "Generated ANIb matrices"
