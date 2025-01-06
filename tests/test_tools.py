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
"""Tests for the pyani_plus/tools.py module.

These tests are intended to be run from the repository root using:

pytest -v

The ``test_bad_*`` tests look at various failures in the command itself,
using a selection of the ``tools.get_*`` functions (which all should behave
the same way).

The ``test_fake_*`` tests call simple bash scripts which echo known strings.
Depending on the tools, this may or may not be parsed as a valid version.

The ``test_find_*`` tests are live and check for the real binary on ``$PATH``.
"""

from pathlib import Path

import pytest

from pyani_plus import tools


def test_bad_path_missing() -> None:
    """Confirm giving an empty path fails."""
    with pytest.raises(
        ValueError, match="Function check_cmd requires a command or full path."
    ):
        tools.get_nucmer("")


def test_bad_path() -> None:
    """Confirm giving an invalid binary path fails."""
    with pytest.raises(RuntimeError, match="/does/not/exist/nucmer is not executable"):
        tools.get_nucmer("/does/not/exist/nucmer")


def test_bad_binary_name() -> None:
    """Confirm giving an invalid binary name fails."""
    with pytest.raises(RuntimeError, match=r"there-is-no-blastn not found on \$PATH"):
        tools.get_blastn("there-is-no-blastn")


def test_non_exec() -> None:
    """Confirm a non-executable binary fails."""
    with pytest.raises(
        RuntimeError,
        match="tests/fixtures/tools/non_executable_script is not executable",
    ):
        tools.get_fastani("tests/fixtures/tools/non_executable_script")


def test_fake_makeblastdb() -> None:
    """Confirm simple makeblastdb version parsing works."""
    info = tools.get_makeblastdb("tests/fixtures/tools/mock_makeblastdb")
    assert info.exe_path == Path("tests/fixtures/tools/mock_makeblastdb").resolve()
    assert info.version == "2.16.0+"

    cmd = Path("tests/fixtures/tools/version_one")  # outputs "version 1.0.0"
    msg = f"Executable exists at {cmd.resolve()} but could not retrieve version"
    with pytest.raises(RuntimeError, match=msg):
        tools.get_makeblastdb(cmd)


def test_fake_blastn() -> None:
    """Confirm simple blastn version parsing works."""
    info = tools.get_blastn("tests/fixtures/tools/mock_blastn")
    assert info.exe_path == Path("tests/fixtures/tools/mock_blastn").resolve()
    assert info.version == "2.16.0+"

    cmd = Path("tests/fixtures/tools/version_one")  # outputs "version 1.0.0"
    msg = f"Executable exists at {cmd.resolve()} but could not retrieve version"
    with pytest.raises(RuntimeError, match=msg):
        tools.get_blastn(cmd)

    cmd = Path("tests/fixtures/tools/just_one")  # outputs "1.0.0"
    msg = f"Executable exists at {cmd.resolve()} but could not retrieve version"
    with pytest.raises(RuntimeError, match=msg):
        tools.get_blastn(cmd)


def test_fake_fastani() -> None:
    """Confirm simple fastani version parsing works."""
    info = tools.get_fastani(
        "tests/fixtures/tools/version_one"
    )  # outputs "version 1.0.0"
    assert info.exe_path == Path("tests/fixtures/tools/version_one").resolve()
    assert info.version == "1.0.0"

    cmd = Path("tests/fixtures/tools/just_one")  # outputs just "1.0.0"
    msg = f"Executable exists at {cmd.resolve()} but could not retrieve version"
    with pytest.raises(RuntimeError, match=msg):
        tools.get_fastani(cmd)


def test_fake_nucmer() -> None:
    """Confirm simple nucmer version parsing works."""
    info = tools.get_nucmer("tests/fixtures/tools/just_one")  # parsed like mummer v4
    assert info.exe_path == Path("tests/fixtures/tools/just_one").resolve()
    assert info.version == "1.0.0"

    info = tools.get_nucmer("tests/fixtures/tools/version_one")  # parsed like mummer v3
    assert info.exe_path == Path("tests/fixtures/tools/version_one").resolve()
    assert info.version == "1.0.0"

    cmd = Path("tests/fixtures/tools/cutting_edge")  # no numerical output
    msg = f"Executable exists at {cmd.resolve()} but could not retrieve version"
    with pytest.raises(RuntimeError, match=msg):
        tools.get_nucmer(cmd)


def test_fake_delta_filter() -> None:
    """Confirm simple delta-filter output parsing works."""
    cmd = Path("tests/fixtures/tools/just_one")  # no numerical output
    msg = f"Executable exists at {cmd.resolve()} but does not seem to be from mummer"
    with pytest.raises(RuntimeError, match=msg):
        tools.get_delta_filter(cmd)


def test_fake_show_coords() -> None:
    """Confirm simple show-coords output parsing works."""
    cmd = Path("tests/fixtures/tools/just_one")  # no numerical output
    msg = f"Executable exists at {cmd.resolve()} but does not seem to be from mummer"
    with pytest.raises(RuntimeError, match=msg):
        tools.get_show_coords(cmd)


def test_fake_show_diff() -> None:
    """Confirm simple show-diff output parsing works."""
    cmd = Path("tests/fixtures/tools/just_one")  # no numerical output
    msg = f"Executable exists at {cmd.resolve()} but does not seem to be from mummer"
    with pytest.raises(RuntimeError, match=msg):
        tools.get_show_diff(cmd)


def test_fake_dnadiff() -> None:
    """Confirm simple dnadiff version parsing works."""
    info = tools.get_nucmer(
        "tests/fixtures/tools/version_one"
    )  # outputs "version 1.0.0"
    assert info.exe_path == Path("tests/fixtures/tools/version_one").resolve()
    assert info.version == "1.0.0"

    info = tools.get_nucmer("tests/fixtures/tools/just_one")  # outputs "1.0.0"
    assert info.exe_path == Path("tests/fixtures/tools/just_one").resolve()
    assert info.version == "1.0.0"

    cmd = Path("tests/fixtures/tools/cutting_edge")  # no numerical output
    msg = f"Executable exists at {cmd.resolve()} but could not retrieve version"
    with pytest.raises(RuntimeError, match=msg):
        tools.get_dnadiff(cmd)


def test_fake_sourmash() -> None:
    """Confirm simple sourmash version parsing works."""
    info = tools.get_sourmash(
        "tests/fixtures/tools/mock_sourmash"
    )  # outputs "sourmash 1.0.0"
    assert info.exe_path == Path("tests/fixtures/tools/mock_sourmash").resolve()
    assert info.version == "1.0.0"

    cmd = Path("tests/fixtures/tools/cutting_edge")  # no numerical output
    msg = f"Executable exists at {cmd.resolve()} but could not retrieve version"
    with pytest.raises(RuntimeError, match=msg):
        tools.get_sourmash(cmd)


def test_find_makeblastdb() -> None:
    """Confirm can find NCBI makeblastdb if on $PATH and determine its version."""
    # At the time of writing this dependency is NOT installed for CI testing
    try:
        info = tools.get_makeblastdb()
    except RuntimeError as err:
        assert str(err) == "makeblastdb not found on $PATH"  # noqa: PT017
    else:
        assert info.exe_path.parts[-1] == "makeblastdb"
        assert info.version.startswith("2.")


def test_find_blastn() -> None:
    """Confirm can find NCBI blastn if on $PATH and determine its version."""
    # At the time of writing this dependency is NOT installed for CI testing
    try:
        info = tools.get_blastn("blastn")
    except RuntimeError as err:
        assert str(err) == "blastn not found on $PATH"  # noqa: PT017
    else:
        assert info.exe_path.parts[-1] == "blastn"
        assert info.version.startswith("2.")


def test_find_fastani() -> None:
    """Confirm can find fastANI on $PATH and determine its version."""
    # At the time of writing this dependency is installed for CI testing
    info = tools.get_fastani()
    assert info.exe_path.parts[-1] == "fastANI"
    assert info.version.startswith("1.")


def test_find_nucmer() -> None:
    """Confirm can find nucmer on $PATH and determine its version."""
    # At the time of writing this dependency is installed for CI testing
    info = tools.get_nucmer()
    assert info.exe_path.parts[-1] == "nucmer"
    assert info.version.startswith("3.")


def test_find_delta_filter() -> None:
    """Confirm can find mummer delta-filter on $PATH."""
    # At the time of writing this dependency is installed for CI testing
    info = tools.get_delta_filter()
    assert info.exe_path.parts[-1] == "delta-filter"
    assert info.version == ""


def test_find_show_coords() -> None:
    """Confirm can find mummer show-coords on $PATH."""
    # At the time of writing this dependency is installed for CI testing
    info = tools.get_show_coords()
    assert info.exe_path.parts[-1] == "show-coords"
    assert info.version == ""


def test_find_show_diff() -> None:
    """Confirm can find mummer show-diff on $PATH."""
    # At the time of writing this dependency is installed for CI testing
    info = tools.get_show_diff()
    assert info.exe_path.parts[-1] == "show-diff"
    assert info.version == ""


def test_find_dnadiff() -> None:
    """Confirm can find dnadiff on $PATH."""
    # At the time of writing this dependency is installed for CI testing
    info = tools.get_dnadiff()
    assert info.exe_path.parts[-1] == "dnadiff"
    assert info.version.startswith("1.")


def test_find_sourmash() -> None:
    """Confirm can find sourmash on $PATH and determine its version."""
    # At the time of writing this dependency is installed for CI testing
    info = tools.get_sourmash()
    assert info.exe_path.parts[-1] == "sourmash"
    assert info.version.startswith("4.")
