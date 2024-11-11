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
"""Code to implement the sourmash Average Nucleotide Identity (ANI) method."""

import platform
import subprocess
import sys
from enum import Enum
from pathlib import Path

import pandas as pd

from pyani_plus import db_orm, tools


class EnumModeSourmash(str, Enum):
    """Enum for the --mode command line argument passed to sourmash."""

    max_containment = "max-containment"  # default
    containment = "containment"


MODE = EnumModeSourmash.max_containment  # constant for CLI default
SCALED = 1000
KMER_SIZE = 31  # default


def parse_compare(compare_file: Path) -> float:
    """Parse sourmash compare .csv output extracting estimated ANI as float.

    :param compare_file: Path to the compare_file.

    Extracts the ANI estimate (which we return in the range 0 to 1).

    Assumes all sourmash comparisons are pairwise, involving one query and
    one reference file. The sourmash compare file should contain a DataFrame
    with two rows and two columns, e.g.:

        genome_1.fna    genome_2.fas
        1.000000        0.997901
        0.997901        1.000000

    To extract the estimated ANI value, we call `df.iloc[0, -1]`
    to obtain the top right corner value.

    Since sourmash *can* produce multi-comparison output when given a list
    of .sig files, it may be necessary to parse the file by extracting
    the upper triangle of values from the sourmash compare output using
    the `np.triu` function.
    """
    compare_results = pd.read_csv(
        Path(compare_file),
        sep=",",
        index_col=False,
    )

    return compare_results.iloc[0, -1]


def compute_pairwise_ani(  # noqa: PLR0913
    uname: platform.uname_result,
    config: db_orm.Configuration,
    query_hash: str,
    query_fasta: Path,
    subject_hash: str,
    subject_fasta: Path,
    cache: Path,
) -> db_orm.Comparison:
    """Run a single sourmash comparison."""
    tool = tools.get_sourmash()
    signature_file = cache / f"{query_hash}_vs_{subject_hash}.sig"
    subprocess.check_call(
        [
            str(tool.exe_path),
            "sketch",
            "dna",
            "-p",
            f"k={config.kmersize},{config.extra}",
            str(query_fasta),
            str(subject_fasta),
            "-o",
            str(signature_file),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if not signature_file.is_file():
        msg = f"ERROR: Missing sourmash signature file {signature_file}"
        sys.exit(msg)
    proc = subprocess.run(
        [
            str(tool.exe_path),
            "compare",
            str(signature_file),
            "--estimate-ani",
            f"--{config.mode}",
            "--quiet",
            "--csv",
            "/dev/stdout",
        ],
        capture_output=True,
        check=True,
        text=True,
    )
    if not proc.stdout:
        msg = f"ERROR: No output from sourmash\n{proc.stderr}"
        sys.exit(msg)

    identity = float(proc.stdout.split("\n")[1].split(",")[1])
    if not (0.0 <= identity <= 1.0):
        msg = f"ERROR: Could not parse sourmash CSV output:\n{proc.stdout}"
        sys.exit(msg)

    return db_orm.Comparison(
        configuration_id=config.configuration_id,
        query_hash=query_hash,
        subject_hash=subject_hash,
        identity=identity,
        uname_system=uname.system,
        uname_release=uname.release,
        uname_machine=uname.machine,
    )
