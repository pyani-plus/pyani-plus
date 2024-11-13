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
from collections.abc import Iterator
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


def prepare_genomes(run: db_orm.Run, cache: Path) -> Iterator[str]:
    """Build the sourmatch sketch signatures in the given directory.

    Yields the FASTA hashes as the databases completed for use with a progress bar.
    """
    config = run.configuration
    if not config.kmersize:
        msg = f"ERROR: sourmash requires a k-mer size, default is {KMER_SIZE}"
        sys.exit(msg)
    if not config.extra:
        msg = f"ERROR: sourmash requires scaled or num, default is scaled={SCALED}"
        sys.exit(msg)
    tool = tools.get_sourmash()
    fasta_dir = Path(run.fasta_directory)
    for entry in run.fasta_hashes:
        fasta_filename = fasta_dir / entry.fasta_filename
        sig_filename = (
            cache / f"{entry.genome_hash}_k={config.kmersize}_{config.extra}.sig"
        )
        if not sig_filename.is_file():
            subprocess.check_call(
                [
                    str(tool.exe_path),
                    "sketch",
                    "dna",
                    "-p",
                    f"k={config.kmersize},{config.extra}",
                    fasta_filename,
                    "-o",
                    sig_filename,
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        yield entry


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
    query_fasta: Path,  # noqa: ARG001
    query_length: int,  # noqa: ARG001
    subject_hash: str,
    subject_fasta: Path,  # noqa: ARG001
    subject_length: int,  # noqa: ARG001
    cache: Path,
) -> db_orm.Comparison:
    """Run a single sourmash comparison."""
    if not config.mode:
        msg = f"ERROR: sourmash requires mode, default is {MODE.value}"
        sys.exit(msg)
    tool = tools.get_sourmash()
    proc = subprocess.run(
        [
            str(tool.exe_path),
            "compare",
            str(cache / f"{query_hash}_k={config.kmersize}_{config.extra}.sig"),
            str(cache / f"{subject_hash}_k={config.kmersize}_{config.extra}.sig"),
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
