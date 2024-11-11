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
"""Code to wrap the fastANI average nucleotide identity method."""

import platform
import subprocess
import sys
from pathlib import Path

from pyani_plus import db_orm, tools

KMER_SIZE = 16  # Default fastANI k-mer size (max 16)
FRAG_LEN = 3000  # Default fastANI fragment length, mapped to our --fragsize
MIN_FRACTION = 0.2  # Default fastANI min fraction
MAX_RATIO_DIFF = 10.0  # Default fastANI maximum ratio difference


def parse_fastani_file(filename: Path) -> tuple[Path, Path, float, int, int]:
    """Parse a single-line fastANI output file extracting key fields as a tuple."""
    with filename.open() as handle:
        line = handle.readline()
    if not line:  # No file content; either run failed or no detectable similarity
        msg = f"Input file {filename} is empty"
        raise ValueError(msg)
    return parse_fastani_data(line)


def parse_fastani_data(text: str) -> tuple[Path, Path, float, int, int]:
    """Parse a single-line fastANI output, extracting key fields as a tuple.

    Return (ref genome, query genome, ANI estimate, orthologous matches,
    sequence fragments) tuple.

    :param filename: Path, path to the input file

    Extracts the ANI estimate (which we return in the range 0 to 1), the
    number of orthologous matches (int), and the number of sequence
    fragments considered from the fastANI output file (int).

    We assume that all fastANI comparisons are pairwise: one query and
    one reference file. The fastANI file should contain a single line.

    fastANI *can* produce multi-line output, if a list of query/reference
    files is given to it.
    """
    line = text.strip().split()
    return (
        Path(line[0]),
        Path(line[1]),
        0.01 * float(line[2]),
        int(line[3]),
        int(line[4]),
    )


def compute_pairwise_ani(  # noqa: PLR0913
    uname: platform.uname_result,
    config: db_orm.Configuration,
    query_hash: str,
    query_fasta: Path,
    subject_hash: str,
    subject_fasta: Path,
    cache: Path,  # noqa: ARG001
) -> db_orm.Comparison:
    """Run a single fastANI comparison."""
    proc = subprocess.run(
        [
            str(tools.get_fastani().exe_path),
            "-q",
            query_fasta,
            "-r",
            subject_fasta,
            "-o",
            "/dev/stdout",
            "--fragLen",
            str(config.fragsize),
            "-k",
            str(config.kmersize),
            "--minFraction",
            str(config.minmatch),
        ],
        capture_output=True,
        check=True,
        text=True,
    )
    if not proc.stdout:
        msg = f"ERROR: No output from fastANI\n{proc.stderr}"
        sys.exit(msg)
    used_query_path, used_subject_path, identity, orthologous_matches, fragments = (
        parse_fastani_data(proc.stdout)
    )
    # Allow for variation in the folder part of the filenames (e.g. relative paths)
    if used_query_path.stem != query_fasta.stem:
        sys.exit(
            f"ERROR: Given --query-fasta {query_fasta} but query in fastANI file contents was {used_query_path}"
        )
    if used_subject_path.stem != subject_fasta.stem:
        sys.exit(
            f"ERROR: Given --subject-fasta {subject_fasta} but subject in fastANI file contents was {used_subject_path}"
        )

    return db_orm.Comparison(
        configuration_id=config.configuration_id,
        query_hash=query_hash,
        subject_hash=subject_hash,
        identity=identity,
        aln_length=round(config.fragsize * orthologous_matches),  # proxy value,
        sim_errors=fragments - orthologous_matches,  # proxy value, not bp,
        cov_query=float(orthologous_matches) / fragments,  # an approximation,
        cov_subject=None,
        uname_system=uname.system,
        uname_release=uname.release,
        uname_machine=uname.machine,
    )
