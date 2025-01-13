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
"""Defines the assorted public command line interface (CLI) arguments we offer.

This is a separate file without any code to allow them to be imported from both
the private and public API definitions without needed to import from each other
(which adds to the start-up time unnecessarily).
"""

from pathlib import Path
from typing import Annotated

import typer

from pyani_plus import FASTA_EXTENSIONS
from pyani_plus.methods.anim import EnumModeANIm
from pyani_plus.workflows import ToolExecutor

# Reused required command line arguments (which have no default)
# --------------------------------------------------------------
# These are named REQ_ARG_TYPE_* short for required-argument type
REQ_ARG_TYPE_DATABASE = Annotated[
    Path,
    typer.Option(
        help="Path to pyANI-plus SQLite3 database",
        show_default=False,
        dir_okay=False,
        file_okay=True,
    ),
]
REQ_ARG_TYPE_OUTDIR = Annotated[
    Path,
    typer.Option(
        help="Output directory",
        show_default=False,
        exists=True,
        dir_okay=True,
        file_okay=False,
    ),
]
REQ_ARG_TYPE_FASTA_DIR = Annotated[
    Path,
    typer.Argument(
        help=f"Directory of FASTA files (extensions {', '.join(sorted(FASTA_EXTENSIONS))})",
        show_default=False,
    ),
]

# Reused optional command line arguments (defined with a default)
# ---------------------------------------------------------------
# These are named OPT_ARG_TYPE_* short for optional-argument type

OPT_ARG_TYPE_RUN_ID = Annotated[
    int | None,
    typer.Option(help="Which run from the database (defaults to latest)"),
]
OPT_ARG_TYPE_RUN_NAME = Annotated[
    # Using None to mean the dynamic default value here:
    str | None,
    typer.Option(
        help="Run name. Default is 'N genomes using METHOD'.", show_default=False
    ),
]
OPT_ARG_TYPE_TEMP_WORKFLOW = Annotated[
    Path | None,
    typer.Option(
        help="Directory to use for temporary workflow coordination files, which"
        " for debugging purposes will not be deleted. For clusters this must be"
        " on a shared drive. Default behaviour is to use a system specified"
        " temporary directory (for the local executor) or a temporary directory"
        " under the present direct (for clusters), and remove this afterwards.",
        rich_help_panel="Debugging",
        show_default=False,
        exists=True,
        dir_okay=True,
        file_okay=False,
    ),
]
OPT_ARG_TYPE_TEMP = Annotated[
    Path | None,
    typer.Option(
        help="Directory to use for intermediate files, which for debugging"
        " purposes will not be deleted. For clusters this must be on a shared"
        " drive. Default behaviour is to use a system specified temporary"
        " directory (specific to the compute-node when using a cluster) and"
        " remove this afterwards.",
        rich_help_panel="Debugging",
        show_default=False,
        exists=True,
        dir_okay=True,
        file_okay=False,
    ),
]
OPT_ARG_TYPE_FRAGSIZE = Annotated[
    int,
    typer.Option(
        help="Comparison method fragment size",
        rich_help_panel="Method parameters",
        min=1,
    ),
]
# fastANI has maximum (and default) k-mer size 16,
# so defined separately with max=16
OPT_ARG_TYPE_KMERSIZE = Annotated[
    int,
    typer.Option(
        help="Comparison method k-mer size",
        rich_help_panel="Method parameters",
        min=1,
    ),
]
OPT_ARG_TYPE_MINMATCH = Annotated[
    float,
    typer.Option(
        help="Comparison method min-match",
        rich_help_panel="Method parameters",
        min=0.0,
        max=1.0,
    ),
]
OPT_ARG_TYPE_ANIM_MODE = Annotated[
    EnumModeANIm,
    typer.Option(
        help="Nucmer mode for ANIm",
        rich_help_panel="Method parameters",
    ),
]
OPT_ARG_TYPE_SOURMASH_SCALED = Annotated[
    int,
    typer.Option(
        help="Sets the compression ratio",
        rich_help_panel="Method parameters",
        min=1,
    ),
]

OPT_ARG_TYPE_CREATE_DB = Annotated[
    # Listing name(s) explicitly to avoid automatic matching --no-create-db
    bool, typer.Option("--create-db", help="Create database if does not exist")
]
OPT_ARG_TYPE_EXECUTOR = Annotated[
    ToolExecutor, typer.Option(help="How should the internal tools be run?")
]
