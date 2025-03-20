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
"""Assorted utility functions used within the pyANI-plus software."""

import gzip
import hashlib
import logging
import os
import shutil
import subprocess
from collections.abc import Iterator
from pathlib import Path
from typing import IO

from pyani_plus import FASTA_EXTENSIONS, log_sys_exit

ASCII_GREATER_THAN = ord(">")  # 64
WHITESPACE = b" \t\r\n"


def fasta_bytes_iterator(
    handle: IO[bytes] | gzip.GzipFile,
) -> Iterator[tuple[bytes, bytes]]:
    """Parse FASTA file in binary (bytes) mode.

    Yields tuples of (description including identifier, sequence).

    >>> with open("tests/fixtures/viral_example/OP073605.fasta", "rb") as handle:
    ...     for title, seq in fasta_bytes_iterator(handle):
    ...         print(title)
    ...         print(f"Length {len(seq)} bp")
    b'OP073605.1 MAG: Bacteriophage sp. isolate 0984_12761, complete genome'
    Length 57793 bp

    Requires a binary mode input handle:

    >>> with open("tests/fixtures/viral_example/OP073605.fasta") as handle:
    ...     for title, seq in fasta_bytes_iterator(handle):
    ...         print(title)
    ...         print(f"Length {len(seq)}bp")
    Traceback (most recent call last):
    ...
    ValueError: Function fasta_bytes_iterator requires a handle in binary mode

    """
    # Follows the logic of the Biopython 1.62 to 1.85 SimpleFastaParser which
    # worked in text (unicode) mode only.
    if handle.read(0) != b"":
        msg = "Function fasta_bytes_iterator requires a handle in binary mode"
        raise ValueError(msg)

    # Skip any text before the first record (e.g. blank lines, comments, header)
    for line in handle:
        if line[0] == ASCII_GREATER_THAN:  # i.e. >
            title = line[1:].rstrip()
            break
    else:
        # no line break, probably an empty file
        return
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines: list[bytes] = []
    for line in handle:
        if line[0] == ASCII_GREATER_THAN:  # i.e. >
            yield title, b"".join(lines).translate(None, WHITESPACE)
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())
    yield title, b"".join(lines).translate(None, WHITESPACE)


def filename_stem(filename: str) -> str:
    """Return the basename (stem) of a filename considering gzip suffix.

    Modified Pathlib stem to consider .gz as well.

    >>> filename_stem("/path/example.fna")
    'example'
    >>> filename_stem("relative/path/example.fna.gz")
    'example'
    """
    if "/" in filename:
        filename = filename.rsplit("/", 1)[1]
    return Path(filename[:-3]).stem if filename.endswith(".gz") else Path(filename).stem


def str_md5sum(text: str, encoding: str = "ascii") -> str:
    """Return the MD5 checksum hash digest of the passed string.

    :param text:  String for hashing

    Behaves like the command line tool ``md5sum``, giving a 32 character
    hexadecimal representation of the MD5 checksum of the given text::

        $ cat tests/fixtures/bacterial_example/NC_002696.fasta.gz | gunzip | md5sum
        f19cb07198a41a4406a22b2f57a6b5e7  -

    In Python:

        >>> with gzip.open(
        ...     "tests/fixtures/bacterial_example/NC_002696.fasta.gz", "rt"
        ... ) as handle:
        ...     text = handle.read()
        >>> str_md5sum(text)
        'f19cb07198a41a4406a22b2f57a6b5e7'

    This particular example would be more consise using the sister function:

        >>> file_md5sum("tests/fixtures/bacterial_example/NC_002696.fasta.gz")
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

    For uncompressed files behaves like the command line tool ``md5sum``,
    giving a 32 character hexadecimal representation of the MD5 checksum
    of the file contents::

        $ md5sum tests/fixtures/viral_example/OP073605.fasta
        5584c7029328dc48d33f95f0a78f7e57  tests/fixtures/viral_example/OP073605.fasta

    In Python:

        >>> file_md5sum("tests/fixtures/viral_example/OP073605.fasta")
        '5584c7029328dc48d33f95f0a78f7e57'

    However, for gzip compressed files, it returns the MD5 of the decompressed
    contents::

        $ cat tests/fixtures/bacterial_example/NC_011916.fas.gz | gunzip | md5sum
        9d72a8fb513cf9cc8cc6605a0ad4e837  -

    In Python:

        >>> file_md5sum("tests/fixtures/bacterial_example/NC_011916.fas.gz")
        '9d72a8fb513cf9cc8cc6605a0ad4e837'

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
        try:
            with gzip.open(fname, "rb") as fhandle:
                for chunk in iter(lambda: fhandle.read(65536), b""):
                    hash_md5.update(chunk)
        except gzip.BadGzipFile:
            with fname.open("rb") as fhandle:
                for chunk in iter(lambda: fhandle.read(65536), b""):
                    hash_md5.update(chunk)
    except FileNotFoundError:
        msg = (
            f"Input {fname} is a broken symlink"
            if fname.is_symlink()
            else f"Input {fname} not found"
        )
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


def check_db(logger: logging.Logger, database: Path | str, create_db: bool) -> None:  # noqa: FBT001
    """Check DB exists, or using create_db=True."""
    msg = f"Checking DB argument '{database}'"
    logger.debug(msg)
    if database != ":memory:" and not create_db and not Path(database).is_file():
        msg = f"Database {database} does not exist, but not using --create-db"
        log_sys_exit(logger, msg)


def check_fasta(logger: logging.Logger, fasta: Path) -> list[Path]:
    """Check fasta is a directory and return list of FASTA files in it."""
    msg = f"Checking FASTA argument '{fasta}'"
    logger.debug(msg)
    if not fasta.is_dir():
        msg = f"FASTA input {fasta} is not a directory"
        log_sys_exit(logger, msg)

    fasta_names: list[Path] = []
    for pattern in FASTA_EXTENSIONS:
        fasta_names.extend(fasta.glob("*" + pattern))
        fasta_names.extend(fasta.glob("*" + pattern + ".gz"))
    if not fasta_names:
        msg = f"No FASTA input genomes under {fasta} with extensions {', '.join(FASTA_EXTENSIONS)}"
        log_sys_exit(logger, msg)

    return fasta_names


def _fmt_cmd(args: list[str]) -> str:
    # Needs to handle spaces with quoting
    return " ".join(f"'{_}'" if " " in str(_) else str(_) for _ in args)


def check_output(logger: logging.Logger, args: list[str]) -> str:
    """Wrap for subprocess.run and log any error.

    Note that if the output is not of interest on success, subprocess.check_call
    would be natural instead. However, the documentation for that recommends not
    to use this with subprocess.PIPE in case the child process blocks.
    """
    code = 0
    output = None
    try:
        return subprocess.check_output(
            [str(_) for _ in args], stderr=subprocess.STDOUT, text=True
        )
    except subprocess.CalledProcessError as err:
        code = err.returncode
        output = err.output

    msg = f"Return code {code} from: {_fmt_cmd(args)}"
    logger.error(msg)  # Is this worth repeating here & on exit?
    logger.error(output)
    if output:
        for line in output.split("\n"):
            if line.upper().startswith("ERROR:"):
                msg += "\n" + line
    # Calling log_sys_exit(X) is equivalent to raising SystemExit(X),
    # where if X is an integer this is the return code, otherwise
    # X is printed to stderr and the return code is 1.
    #
    # i.e. We cannot raise SystemExit with a custom message and code
    # This makes testing with pytest a little harder! Therefore opting
    # to raise a message which we can test, and settle for return code 1.
    log_sys_exit(logger, msg)
    # mypy wants to see a return:
    return ""  # pragma: nocover


def stage_file(
    logger: logging.Logger,
    input_filename: Path,
    staged_filename: Path,
    *,
    decompress: bool = True,
) -> None:
    """Prepare a symlink or decompressed copy of the given file.

    This is used to avoid spaces and other problematic characters in input
    FASTA filenames by staging a symlink named ``<md5>.fasta``, and to handle
    gzipped input transparently by decompressing to ``<md5>.fasta`` instead.

    We use this helper function for calling ``mummer`` and ``makeblastdb``, but
    it is not needed for ``fastANI`` or ``sourmash``.
    """
    if not input_filename.is_file():
        msg = f"Missing input file {input_filename}"
        log_sys_exit(logger, msg)
    if staged_filename.is_file():
        # This could be a race condition?
        # Perhaps if resume with explicit temp directory given?
        msg = f"Intermediate file {staged_filename} already exists!"
        log_sys_exit(logger, msg)

    if input_filename.suffix == ".gz" and decompress:
        # Decompress to the given temporary filename
        with (
            gzip.open(input_filename, "rb") as f_in,
            staged_filename.open("wb") as f_out,
        ):
            shutil.copyfileobj(f_in, f_out)
    else:
        # Make a symlink pointing to the original
        staged_filename.symlink_to(input_filename)
