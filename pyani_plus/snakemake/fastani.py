"""Calling snakemake scheduler from Python to generate
target outputs (delta files) for aniM analysis.
"""

# Set Up (importing libraries)
from importlib import resources as impresources
from itertools import permutations
from pathlib import Path

from snakemake.api import SnakemakeApi, _get_executor_plugin_registry
from snakemake.settings import ConfigSettings, DAGSettings, ResourceSettings

from pyani_plus import workflows

# config_args = {
#     "outdir": "../../issue_1/output5",
#     "indir": "../../issue_1/input",
#     "cores": 8,
#     "fragLen": 3000,
#     "kmerSize": 16,
#     "minFrac": 0.2,
# }


# targets = [
#     "../../issue_1/output5/NC_002696_vs_NC_010338.fastani",
#     "../../issue_1/output5/NC_002696_vs_NC_011916.fastani",
#     "../../issue_1/output5/NC_002696_vs_NC_014100.fastani",
#     "../../issue_1/output5/NC_010338_vs_NC_002696.fastani",
#     "../../issue_1/output5/NC_010338_vs_NC_011916.fastani",
#     "../../issue_1/output5/NC_010338_vs_NC_014100.fastani",
# ]


def run_workflow(targetset, config_args):
    """Runs the snakemake_fastani workflow for the passed pairwise comparisons"""
    targetset = [str(_) for _ in targetset]
    print(f"{targetset=}")

    # Path to anim snakemake file
    # snakefile = Path("../workflows/snakemake_fastani.smk")
    snakefile = impresources.files(workflows) / "snakemake_fastani.smk"

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


# run_workflow(targets,config_args)
