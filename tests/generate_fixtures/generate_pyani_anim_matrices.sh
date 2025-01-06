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

echo "This is intended to be used with pyANI v0.3, we have:"
pyani --version
# Abort if not matched (via set -e above)
pyani --version 2>&1 | grep "pyani version: 0\.3\." > /dev/null

# Remove pre-existing fixtures before regenerating new ones.
# This is to help with if and when we change the
# example sequences being used.
echo "Removing pre-existing fixtures..."
for file in ../fixtures/viral_example/matrices/ANIm_*.tsv; do
    rm "$file"
done

# Make a temp subdir
rm -rf tmp_pyani_anim
mkdir tmp_pyani_anim
cd tmp_pyani_anim

echo "Creating pyANI database"
pyani createdb

echo "Setting up input genomes"
for FASTA in ../../../tests/fixtures/viral_example/*.f*; do
    ln -s "$FASTA" "${FASTA##*/}"
done

echo "Indexing with pyani..."
#index with pyANI v0.3
pyani index -i .

#Update the label file so that the columns in the matrices correspond to the stem
awk -F'\t' 'OFS="\t" {$3 = $2; print}' labels.txt > labels_tmp.txt && mv labels_tmp.txt labels.txt

echo "Run ANIm comparisions..."
pyani anim -i . -o output -v -l output/log.log --name "generate fixtures" --labels labels.txt --classes classes.txt

echo "Generate matrices..."
pyani report -v -o . --formats=stdout --run_matrices 1

echo "Collecting output for test fixtures..."
for file in matrix_*_1.tab; do
    # Renaming matrix_XXX_1.tab to ANIb_XXX.tsv
    name=${file%_1.tab}
    name=${name#matrix_}
    # Removing :n suffices on MD5 hash labels
    sed "s/:[0-9]//g" "$file" > "../../fixtures/viral_example/matrices/ANIm_${name}.tsv"
done

#Remove temp subdir
cd ..
rm -rf tmp_pyani_anim
echo "Generated ANIm matrices"
