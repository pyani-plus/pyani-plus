#!/bin/bash
set -euo pipefail
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

echo "This is intended to be used with pyani-plust v1+, we have:"
pyani-plus --version
# Abort if not matched (via set -e above)
pyani-plus --version 2>&1 | grep "pyANI-plus 1\." > /dev/null

# Remove pre-existing fixtures before regenerating new ones.
# This is to help with if and when we change the
# example sequences being used.
echo "Removing pre-existing fixtures..."
rm ../fixtures/viral_example/matrices/lzani_*.tsv || :

# Make a temp subdir
rm -rf tmp_pyani_lzani
mkdir tmp_pyani_lzani

echo "Run lzani comparisions..."
pyani-plus lzani --create-db -d tmp_pyani_lzani.db ../fixtures/viral_example/

echo "Generate matrices..."
pyani-plus export-run -d tmp_pyani_lzani.db -o tmp_pyani_lzani -r 1

echo "Collecting output for test fixtures..."
cp tmp_pyani_lzani/*.tsv ../fixtures/viral_example/matrices/

#Remove temp subdir and database
rm -rf tmp_pyani_lzani
rm tmp_pyani_lzani.db
echo "Generated lzani matrices"
