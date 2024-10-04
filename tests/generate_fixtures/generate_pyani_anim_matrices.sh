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

# Make a temp subdir
mkdir tmp_pyani
cd tmp_pyani

# Create pyanidb
pyani createdb

# Copy genomes to tmp_pyani subdir (macOS)
for FASTA in ../../../tests/fixtures/viral_example/*.f*; do
    ln -s $FASTA ${FASTA##*/};
done

#index with pyANI v0.3
pyani index -i .

#Update the label file so that the columns in the matrices correspond to the genome md5 hash
awk -F'\t' 'OFS="\t" {$3 = $1; print}' labels.txt > labels_tmp.txt && mv labels_tmp.txt labels.txt

#Run comparisions
pyani anim -i . -o output -v -l output/log.log --name "genearte fixtures" --labels labels.txt --classes classes.txt

#Generate matrices
pyani report -v -o ../../../tests/fixtures/anim/matrices/ --formats=stdout --run_matrices 1

#Change extensions files
for file in ../../../tests/fixtures/anim/matrices/*; do
    mv "$file" "${file%.tab}.tsv"
done

#Remove temp subdir
cd ..
rm -rf tmp_pyani
