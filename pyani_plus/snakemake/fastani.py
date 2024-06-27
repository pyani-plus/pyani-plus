"""Interface to snakemake that generates target outputs (delta files) for fastANI analysis."""

# Set Up (importing libraries)
from importlib import resources as impresources

from snakemake.api import (  # type: ignore  # noqa: PGH003
    ConfigSettings,
    DAGSettings,
    ResourceSettings,
    SnakemakeApi,
)

from pyani_plus import workflows


def run_workflow(targetset: list[str], config_args: dict) -> None:
    """Run snakemake_fastani workflow for passed pairwise comparisons."""
    targetset = [str(_) for _ in targetset]

    # Path to anim snakemake file
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
            ),
        )
        dag_api.execute_workflow()
