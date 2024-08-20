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


def test_fake_blastn() -> None:
    """Confirm simple blastn version parsing works."""
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
    """Confirm simple fastani version parsing works."""
    info = tools.get_nucmer("tests/fixtures/tools/just_one")  # parsed like mummer v4
    assert info.exe_path == Path("tests/fixtures/tools/just_one").resolve()
    assert info.version == "1.0.0"

    info = tools.get_nucmer("tests/fixtures/tools/version_one")  # parsed like mummer v3
    assert info.exe_path == Path("tests/fixtures/tools/version_one").resolve()
    assert info.version == "1.0.0"


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
