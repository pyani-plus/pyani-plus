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
"""Assorted utility functions used within the pyANI-plus software."""

import hashlib
import os
import sys
from pathlib import Path

from pyani_plus import FASTA_EXTENSIONS


def str_md5sum(text: str, encoding: str = "ascii") -> str:
    """Return the MD5 checksum hash digest of the passed string.

    :param text:  String for hashing

    Behaves like the command line tool ``md5sum``, giving a 32 character
    hexadecimal representation of the MD5 checksum of the file contents::

        $ md5sum tests/fixtures/sequences/NC_002696.fasta
        f19cb07198a41a4406a22b2f57a6b5e7  tests/fixtures/sequences/NC_002696.fasta

    In Python:

        >>> with open("tests/fixtures/sequences/NC_002696.fasta") as handle:
        ...     text = handle.read()
        >>> str_md5sum(text)
        'f19cb07198a41a4406a22b2f57a6b5e7'

    This particular example would be more consise using the sister function:

        >>> file_md5sum("tests/fixtures/sequences/NC_002696.fasta")
        'f19cb07198a41a4406a22b2f57a6b5e7'

    The MD5 checksum is used in pyANI-plus on input FASTA format sequence files.
    This helper function is to avoid repetitive code and associated warnings
    from code checking tools, and is used in our test suite.
    """
    # We're ignoring the linter warning as not using MD5 for security:
    # S324 Probable use of insecure hash functions in `hashlib`: `md5`
    return hashlib.md5(text.encode(encoding)).hexdigest()  # noqa: S324


def file_md5sum(filename: Path | str) -> str:
    """Return the MD5 checksum hash digest of the passed file contents.

    :param filename:  Path or string, path to file for hashing

    Behaves like the command line tool ``md5sum``, giving a 32 character
    hexadecimal representation of the MD5 checksum of the file contents::

        $ md5sum tests/fixtures/viral_example/OP073605.fasta
        5584c7029328dc48d33f95f0a78f7e57  tests/fixtures/viral_example/OP073605.fasta

    In Python:

        >>> file_md5sum("tests/fixtures/viral_example/OP073605.fasta")
        '5584c7029328dc48d33f95f0a78f7e57'

    This is used in pyANI-plus on input FASTA format sequence files, to give a
    fingerprint of the file contents allowing us to cache and reused comparison
    results even when the sequence files are renamed or moved. Note any change
    to the file contents (e.g. editing a description) will change the checksum.
    """
    fname = Path(filename)  # ensure we have a Path object
    # We're ignoring the linter warning as not using MD5 for security:
    # S324 Probable use of insecure hash functions in `hashlib`: `md5`
    hash_md5 = hashlib.md5()  # noqa: S324
    try:
        with fname.open("rb") as fhandle:
            for chunk in iter(lambda: fhandle.read(65536), b""):
                hash_md5.update(chunk)
    except FileNotFoundError:
        msg = f"Input file {fname} is not a file or symlink"
        raise ValueError(msg) from None

    return hash_md5.hexdigest()


def available_cores() -> int:
    """How many CPU cores/threads are available to use."""
    try:
        # This will take into account SLURM limits,
        # so don't need to check $SLURM_CPUS_PER_TASK explicitly.
        # Probably don't need to check $NSLOTS on SGE either.
        available = len(os.sched_getaffinity(0))  # type: ignore[attr-defined]
    except AttributeError:
        # Unavailable on macOS or Windows, use this instead
        # Can return None (but under what circumstances?)
        cpus = os.cpu_count()
        if not cpus:
            msg = "Cannot determine CPU count"  # pragma: no cover
            raise RuntimeError(msg) from None  # pragma: no cover
        available = cpus
    return available


def check_db(database: Path | str, create_db: bool) -> None:  # noqa: FBT001
    """Check DB exists, or using create_db=True."""
    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"ERROR: Database {database} does not exist, but not using --create-db"
        sys.exit(msg)


def check_fasta(fasta: Path) -> list[Path]:
    """Check fasta is a directory and return list of FASTA files in it."""
    if not fasta.is_dir():
        msg = f"ERROR: FASTA input {fasta} is not a directory"
        sys.exit(msg)

    fasta_names: list[Path] = []
    for pattern in FASTA_EXTENSIONS:
        fasta_names.extend(fasta.glob("*" + pattern))
    if not fasta_names:
        msg = f"ERROR: No FASTA input genomes under {fasta} with extensions {', '.join(FASTA_EXTENSIONS)}"
        sys.exit(msg)

    return fasta_names
