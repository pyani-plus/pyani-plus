#!/usr/bin/env bash

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

#create db
pyani createdb --force

#index with pyANI v0.3
pyani index -i ../../tests/fixtures/viral_example/

#Update the label file so that the columns in the matrices correspond to the genome stem files
awk -F'\t' 'OFS="\t" {$3 = $2; print}' ../../tests/fixtures/viral_example/labels.txt > ../../tests/fixtures/viral_example/labels_tmp.txt && mv ../../tests/fixtures/viral_example/labels_tmp.txt ../../tests/fixtures/viral_example/labels.txt

#Run comparisions
pyani anim -i ../../tests/fixtures/viral_example/ -o ../../tests/fixtures/viral_example/output -v -l ../../tests/fixtures/viral_example/output/log.log --name "genearte fixtures" --labels ../../tests/fixtures/viral_example/labels.txt --classes ../../tests/fixtures/viral_example/classes.txt

#Generate matrices
pyani report -v --runs -o ../../tests/fixtures/anim/matrices/ --formats=stdout --run_matrices 1

#Change extensions files
for file in ../../tests/fixtures/anim/matrices/*; do
    mv "$file" "${file%.tab}.tsv"
done

#Remove unwanted files (eg. pyANI output)
rm -r ../../tests/fixtures/viral_example/output
rm -r ../../tests/fixtures/viral_example/*.md5
rm -r ../../tests/fixtures/viral_example/*.txt
