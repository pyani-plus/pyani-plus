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

import logging
from pathlib import Path

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
GRAPHICS_FORMATS = ("tsv", "png", "jpg", "svgz", "pdf")  # note no dots!
PROGRESS_BAR_COLUMNS = [
    TextColumn("[progress.description]{task.description}"),
    BarColumn(),
    TaskProgressColumn(),
    # Removing TimeRemainingColumn() from defaults, replacing with:
    TimeElapsedColumn(),
    # Add this last as have some out of N and some out of N^2:
    MofNCompleteColumn(),
]


def setup_logger(log_folder: Path | None, name: str) -> logging.Logger:
    """Return a file-based logger alongside a Rich console logger.

    The file logger defaults to DEBUG level, but the less verbose INFO for the terminal.
    """
    logger = logging.getLogger(f"{__package__}.{name}")
    logger.setLevel(logging.DEBUG)
    if logger.hasHandlers():  # remove all previous handlers to avoid duplicate entries
        logger.handlers.clear()

    if log_folder == Path("-"):
        logger.debug("Not logging to file.")
        return logger

    filename = (log_folder / (name + ".log")) if log_folder else Path(name + ".log")
    file_handler = logging.FileHandler(filename, mode="a")
    file_handler.setLevel(logging.DEBUG)

    fmt = "%(asctime)s %(levelname)9s %(filename)21s:%(lineno)-3s | %(message)s"
    formatter = logging.Formatter(fmt=fmt, datefmt="%Y-%m-%d %H:%M:%S")
    file_handler.setFormatter(formatter)

    logger.addHandler(file_handler)

    msg = f"Logging to {filename}"
    logger.info(msg)  # Want this to appear on the terminal

    return logger
