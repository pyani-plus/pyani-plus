# Set Up (importing libraries)
from importlib import resources as impresources
from itertools import permutations
from pathlib import Path

from snakemake.api import SnakemakeApi, ConfigSettings, DAGSettings, ResourceSettings

from pyani_plus import workflows


def check_input_stems(indir):
    """Check input files agans approved list of extenions.
    If duplicate stems with approved extenions are present
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
        raise ValueError(
            f"Duplicated stems found for {list(set(duplicates))}. Please investigate."
        )

    return input_files


class SnakemakeRunner:
    def __init__(self, workflow_filename):
        self.workflow_filename = workflow_filename

    def run_workflow(self, targetset, config_args):
        """Runs the snakemake workflow"""
        targetset = [str(_) for _ in targetset]
        print(f"{targetset=}")

        # Path to anim snakemake file
        snakefile = impresources.files(workflows) / self.workflow_filename

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
            dag_api.unlock()
            dag_api.execute_workflow()
