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

"""pyANI-plus.

This is ``pyANI-plus``, an application and Python module for whole-genome
classification of microbes using Average Nucleotide Identity (ANI) and similar
methods. It is a reimplemented version of ``pyani`` with support for
additional schedulers and methods.
"""

from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)

__version__ = "0.0.1"

# The following are assorted centrally defined constants:
FASTA_EXTENSIONS = {".fasta", ".fas", ".fna"}  # we'll consider .fasta.gz etc too
GRAPHICS_FORMATS = ("tsv", "png", "jpg", "svg", "pdf")  # note no dots!
PROGRESS_BAR_COLUMNS = [
    TextColumn("[progress.description]{task.description}"),
    BarColumn(),
    TaskProgressColumn(),
    # Removing TimeRemainingColumn() from defaults, replacing with:
    TimeElapsedColumn(),
    # Add this last as have some out of N and some out of N^2:
    MofNCompleteColumn(),
]
