"""Calling snakemake scheduler from Python to generate
target outputs (delta files) for dnadiff analysis.
"""

# Set Up (importing libraries)
from importlib import resources as impresources
from itertools import permutations
from pathlib import Path

from snakemake.api import SnakemakeApi, _get_executor_plugin_registry
from snakemake.api import ConfigSettings, DAGSettings, ResourceSettings

from pyani_plus import workflows


def run_workflow(targetset, config_args):
    """Runs the snakemake_fastani workflow for the passed pairwise comparisons"""
    targetset = [str(_) for _ in targetset]
    print(f"{targetset=}")

    # Path to anim snakemake file
    snakefile = impresources.files(workflows) / "snakemake_dnadiff.smk"

    # Use the defined workflow from the Python API
    with SnakemakeApi() as snakemake_api:
        workflow_api = snakemake_api.workflow(
            snakefile=snakefile,
            resource_settings=ResourceSettings(cores=config_args["cores"]),
            config_settings=ConfigSettings(
                config=config_args,
            ),
        )
        dag_api = workflow_api.dag(
            dag_settings=DAGSettings(
                targets=targetset,
            )
        )
        dag_api.execute_workflow()

    pass
