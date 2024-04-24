"""Calling snakemake scheduler from Python to generate 
target outputs (delta files) for aniM analysis. 
"""


# Set Up (importing libraries)
from pathlib import Path
from snakemake.api import SnakemakeApi, _get_executor_plugin_registry
from snakemake.settings import ConfigSettings, DAGSettings, ResourceSettings
from itertools import permutations


# Define a subset of target files to generate
def get_target_files(input_genomes):

    """Return list of target files. 

    The target files include delta-filter files, for both
    forward and reverse pairwise comparisions. 

    :param input_genomes: Path to input genomes files
    """

    #Path to input files
    datadir = Path(input_genomes)
    filenames = sorted(datadir.glob("*"))

    #Just incase there are other hidden files in the input directory...
    #We can construct target files for given FASTA files
    FASTA_extensions = ['.fna', '.fasta']



    comparisions = list(permutations([fname.stem for fname in filenames if fname.suffix in FASTA_extensions], 2))
    target_file = [_[0] + '_vs_' + _[1] + '.filter' for _ in comparisions]


    return target_file




if __name__ == "__main__":
    # Set Up (importing libraries)
    from snakemake.api import SnakemakeApi
    from snakemake.settings import ConfigSettings, DAGSettings, ResourceSettings

    # Path to anim snakemake file
    snakefile = Path("snakemake_workflows/snakemake_anim.smk")

    # Define arguments to pass to the snakefile (eg. input, output, cores, mode)
    config_args = {
        "outdir": "../issue_1/output",
        "indir": "../issue_1/input",
        "cores": 8,
        "mode": 'maxmatch'
    }

    target_files = get_target_files(config_args['indir'])

    # Use the defined workflow from the Python API
    with SnakemakeApi() as snakemake_api:
        workflow_api = snakemake_api.workflow(snakefile=snakefile,
                                              resource_settings=ResourceSettings(cores=config_args["cores"]),
                                              config_settings=ConfigSettings(
                                                  config=config_args,
                                              )
                                              )
        dag_api = workflow_api.dag(
            dag_settings = DAGSettings(
                targets=[config_args['outdir'] + '/' + _ for _ in target_files],  
            )
        )
        dag_api.execute_workflow()