"""Interface to snakemake to generate target outputs (delta files) for ANIm."""

# Set Up (importing libraries)
from importlib import resources as impresources

# Although `from snakemake.settings import ConfigSettings` works fine for <8.14
# the new import choice is required after that point but is backwards-compatible.
from snakemake.api import (  # type: ignore # noqa:PGH003
    ConfigSettings,
    DAGSettings,
    ResourceSettings,
    SnakemakeApi,
)

from pyani_plus import workflows


# FOR DEMONSTRATION PURPOSES ONLY - UNTESTED AND MAY NOT WORK
def run_workflow(targetset: list[str], config_args: dict) -> None:
    """Run snakemake_anim workflow for the passed pairwise comparisons.

    - targetset: list of target output files to trigger snakemake rules
    - config_args:

    This function should be able to take a collection of targets of arbitrary
    size, and pass them to the workflow for distribution, e.g.

    run_workflow(["A_vs_B.filter"])

    mytargets = (["X_vs_Y.filter", "Y_vs_X.filter", "G_vs_M.filter"])
    run_workflow(mytargets)

    We should expect this to be the standard way to run the workflow - passing
    a list of targets to the final rule in the dependency graph - as we will
    be using it to pick up missed comparisons/compare only unanalysed genomes,
    etc.
    """
    targetset = [str(_) for _ in targetset]

    # Path to anim snakemake file
    snakefile = impresources.files(workflows) / "snakemake_anim.smk"

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
