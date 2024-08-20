# (c) The University of Strathclyde 2024-present
# Author: Peter Cock
#
# Contact:
# peter.cock@strath.ac.uk
#
# Peter Cock,
# Strathclyde Institute of Pharmaceutical and Biomedical Sciences
# The University of Strathclyde
# 161 Cathedral Street
# Glasgow
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# (c) The University of Strathclyde 2024-present
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
"""Assorted functions for capturing command line tool version information.

The pyANI-plus software relies on multiple separate command line tools
such as ``nucmer`` from mummer3 to computer average nucleotide identity.
This module contains helper functions to get the version of such tools,
typically by executing them and parsing the stdout. These can be used for
logging, as the exact versions used can sometimes be important for full
reproducibility of results.
"""

import re
import shutil
import subprocess
from pathlib import Path
from typing import NamedTuple


class ExternalToolData(NamedTuple):
    """Convenience struct for tool path and version information."""

    exe_path: Path
    version: str


def check_cmd(cmd: str | Path) -> Path:
    """Return full Path of given command, or raise RuntimeError."""
    if not cmd:
        msg = "Function check_cmd requires a command or full path."
        raise ValueError(msg)

    exe_str = shutil.which(cmd)
    if not exe_str:
        msg = (
            f"{cmd} is not executable"
            if "/" in str(cmd)
            else f"{cmd} not found on $PATH"
        )
        raise RuntimeError(msg)

    # Turn into an absolute path in case later use a temp working dir:
    return Path(exe_str).absolute()


def _get_path_and_version_output(
    cmd: str | Path, args: list[str] | None = None
) -> tuple[Path, str]:
    """Determine path of command, run it with args, capture combined stdout and stderr."""
    # Might later need to add check=True as an optional argument,
    # e.g. NCBI legacy blast doesn't have a version option and uses return code 1.
    exe_path = check_cmd(cmd)
    result = subprocess.run(
        [str(exe_path), *(args if args else [])],
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=True,
        text=True,
    )
    assert isinstance(result.stdout, str)  # noqa: S101
    return exe_path, result.stdout


def get_blastn(cmd: str | Path = "blastn") -> ExternalToolData:
    """Return NCBI BLAST+ blastn path and version as a named tuple.

    We expect the tool to behave as follows:

    .. code-block:: bash

        $ blastn -version
        blastn: 2.16.0+
         Package: blast 2.16.0, build Aug  6 2024 15:58:44

    Here the function would return the binary path and "2.16.0+" as the version.

    Raises a RuntimeError if the command cannot be found, run, or if the
    version cannot be inferred - this likely indicates a dramatically
    different version of the tool with different behaviour.
    """
    exe_path, output = _get_path_and_version_output(cmd, ["-version"])

    match = re.search(r"(?<=blastn:\s)[0-9\.]*\+", output)
    version = match.group().strip() if match else None
    if not version:
        msg = f"Executable exists at {exe_path} but could not retrieve version"
        raise RuntimeError(msg)

    return ExternalToolData(exe_path, version)


def get_fastani(cmd: str | Path = "fastANI") -> ExternalToolData:
    """Return FastANI path and version as a named tuple.

    We expect the tool to behave as follows:

    .. code-block:: bash

        $ fastANI -v
        version 1.33

    Here the function would return the binary path and "1.33" as the version.

    Raises a RuntimeError if the command cannot be found, run, or if the
    version cannot be inferred - this likely indicates a dramatically
    different version of the tool with different behaviour.
    """
    exe_path, output = _get_path_and_version_output(cmd, ["-v"])

    match = re.search(r"(?<=version\s)[0-9\.]*", output)
    version = match.group().strip() if match else None
    if not version:
        msg = f"Executable exists at {exe_path} but could not retrieve version"
        raise RuntimeError(msg)

    return ExternalToolData(exe_path, version)


def get_nucmer(cmd: str | Path = "nucmer") -> ExternalToolData:
    """Return NUCmer path and version as a named tuple.

    :param cmd:  path to NUCmer executable.

    We expect NUCmer v3 to return a string on STDERR as

    .. code-block:: bash

        $ nucmer -V
        nucmer
        NUCmer (NUCleotide MUMmer) version 3.1

    In this case the return value would be the full path of the
    command, and "3.1" as the version.

    For v4, the output is shorter and to STDOUT:

    .. code-block:: bash

        $ nucmer -V
        4.0.0rc1

    Here the function would return "4.0.0.rc1" as the version.

    Raises a RuntimeError if the command cannot be found, run, or if the
    version cannot be inferred - this likely indicates a dramatically
    different version of the tool with different behaviour.
    """
    exe_path, output = _get_path_and_version_output(cmd, ["-V"])

    version = None
    # Try nucmer v3 style
    match = re.search(r"(?<=version\s)[0-9\.]*", output)
    version = match.group().strip() if match else None
    if not version:
        # Try nucmer v4 style
        match = re.search(r"[0-9a-z\.]*", output)
        version = match.group().strip() if match else None
    if not version:
        msg = f"Executable exists at {exe_path} but could not retrieve version"
        raise RuntimeError(msg)

    return ExternalToolData(exe_path, version)