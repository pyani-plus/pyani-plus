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
"""Snakemake scheduler for pairwise comparisons for ANI analysis."""

import importlib.resources as impresources
from pathlib import Path

from snakemake.api import ConfigSettings, DAGSettings, ResourceSettings, SnakemakeApi

from pyani_plus import workflows


def check_input_stems(indir: str) -> dict[str, Path]:
    """Check input files against approved list of extensions.

    If duplicate stems with approved extensions are present
    raise a ValueError.
    """
    extensions = [".fasta", ".fas", ".fna"]
    stems = [_.stem for _ in Path(indir).glob("*") if _.suffix in extensions]

    if len(stems) == len(set(stems)):
        input_files = {
            _.stem: _ for _ in Path(indir).glob("*") if _.suffix in extensions
        }
    else:
        duplicates = [
            item for item in stems if stems.count(item) > 1 and item in set(stems)
        ]
        msg = (
            f"Duplicated stems found for {sorted(set(duplicates))}. Please investigate."
        )
        raise ValueError(msg)

    return input_files


class SnakemakeRunner:
    """Execute Snakemake workflows for pairwise comparisons."""

    def __init__(self, workflow_filename: str) -> None:
        """Initialise SnakemakeRunner with SnakeMake workflow filename."""
        self.workflow_filename = workflow_filename

    def run_workflow(
        self, targetset: list, config_args: dict, workdir: Path | None = None
    ) -> None:
        """Run the snakemake workflow.

        :param targets: (list) List of target files that the workflow should return.
        :config_args: (dict) Dictionary of configuration arguments to pass to the smk workflow.
        """
        targetset = [str(_) for _ in targetset]

        # Path to anim snakemake file
        snakefile = impresources.files(workflows) / self.workflow_filename

        # Use the defined workflow from the Python API
        with SnakemakeApi() as snakemake_api:
            workflow_api = snakemake_api.workflow(
                snakefile=snakefile,
                resource_settings=ResourceSettings(cores=config_args["cores"]),
                config_settings=ConfigSettings(config=config_args),
                workdir=workdir,
            )
            dag_api = workflow_api.dag(
                dag_settings=DAGSettings(targets=targetset),
            )
            dag_api.unlock()
            dag_api.execute_workflow()
