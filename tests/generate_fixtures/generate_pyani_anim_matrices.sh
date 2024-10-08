#!/bin/bash
set -euo pipefail
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

echo "This is intended to be used with pyANI version 0.2.13.1, we have:"
average_nucleotide_identity.py --version
# Abort if not matched (via set -e above)
average_nucleotide_identity.py --version | grep "pyani 0\.2\." > /dev/null

# Make a temp subdir
mkdir tmp_pyani_ANIm
cd tmp_pyani_ANIm

echo "Setting up input genomes"
for FASTA in ../../../tests/fixtures/viral_example/*.f*; do
    ln -s "$FASTA" "${FASTA##*/}"
done

# This gives use MD5 spaces filename with extension:
#     md5sum -- *.f*
# This drops the extension (assuming only one dot present):
#     cut -f 1 -d "."
# This swaps the field order and turns it into a sed search-and-replace command:
#     awk '{print "s/" $2 "/" $1 "/g"}'
md5sum -- *.f* | cut -f 1 -d "." | awk '{print "s/" $2 "/" $1 "/g"}' > rename.txt

echo "Run ANIm comparisions..."
mkdir output
average_nucleotide_identity.py -m ANIm -i . -o output -v -l ANIm.log --force # --labels labels.txt --classes classes.txt

echo "Collecting output for test fixtures..."
# The output names here are based on pyANI v0.3 conventions
SORTING="import sys; import pandas as pd; pd.read_csv(sys.stdin, sep='\t', index_col=0).sort_index(axis=0).sort_index(axis=1).to_csv(sys.stdout, sep='\t')"
sed -f rename.txt output/ANIm_percentage_identity.tab | python -c "$SORTING" > ../../fixtures/ANIm/matrices/matrix_identity.tsv
sed -f rename.txt output/ANIm_alignment_coverage.tab | python -c "$SORTING" > ../../fixtures/ANIm/matrices/matrix_coverage.tsv
sed -f rename.txt output/ANIm_alignment_lengths.tab | python -c "$SORTING" > ../../fixtures/ANIm/matrices/matrix_aln_lengths.tsv
sed -f rename.txt output/ANIm_hadamard.tab | python -c "$SORTING" > ../../fixtures/ANIm/matrices/matrix_hadamard.tsv
sed -f rename.txt output/ANIm_similarity_errors.tab | python -c "$SORTING" > ../../fixtures/ANIm/matrices/matrix_sim_errors.tsv

#Remove temp subdir
cd ..
rm -rf tmp_pyani_ANIm
echo "Generated ANIm matrices"
