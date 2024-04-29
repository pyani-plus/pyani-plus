"""Calling snakemake scheduler from Python to generate 
target outputs (delta files) for aniM analysis. 
"""

# Set Up (importing libraries)
from pathlib import Path
from snakemake.api import SnakemakeApi, _get_executor_plugin_registry
from snakemake.settings import ConfigSettings, DAGSettings, ResourceSettings
from itertools import permutations



config_args = {
        "outdir": "../issue_1/output3",
        "indir": "../issue_1/input",
        "cores": 8,
        "mode": 'maxmatch'
    }


targets = ['../issue_1/output3/NC_002696_vs_NC_010338.filter', 
           '../issue_1/output3/NC_002696_vs_NC_011916.filter', 
           '../issue_1/output3/NC_002696_vs_NC_014100.filter', 
           '../issue_1/output3/NC_010338_vs_NC_002696.filter', 
           '../issue_1/output3/NC_010338_vs_NC_011916.filter', 
           '../issue_1/output3/NC_010338_vs_NC_014100.filter']

# FOR DEMONSTRATION PURPOSES ONLY - UNTESTED AND MAY NOT WORK
def run_workflow(targetset):
    """Runs the snakemake_anim workflow for the passed pairwise comparisons

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

    # Path to anim snakemake file
    snakefile = Path("snakemake_workflows/snakemake_anim.smk")


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
                targets=targetset,  
            )
        )
        dag_api.execute_workflow()

    pass

run_workflow(targets)


config_args = {
        "outdir": "../issue_1/output4",
        "indir": "../issue_1/input",
        "cores": 8,
        "mode": 'maxmatch'
    }


# FOR DEMONSTRATION PURPOSES ONLY - UNTESTED AND MAY NOT WORK
def run_workflow_dir(dirpath):
    """Runs the snakemake_anim workflow on all files in the passed directory

  

    We don't necessarily need this for running pyani-plus, but it may be useful.
    Mostly we will expect to have to filter out some pre-calculated comparisons,
    in normal use.
    """

    # Path to anim snakemake file
    snakefile = Path("snakemake_workflows/snakemake_anim.smk")

    #Path to input files
    datadir = Path(dirpath)
    filenames = sorted(datadir.glob("*"))

    #Just incase there are other hidden files in the input directory...
    #We can construct target files for given FASTA files
    FASTA_extensions = ['.fna', '.fasta']


    comparisions = list(permutations([fname.stem for fname in filenames if fname.suffix in FASTA_extensions], 2))
    targetset = [config_args['outdir'] + '/' + _[0] + '_vs_' + _[1] + '.filter' for _ in comparisions]


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
                targets=targetset,  
            )
        )
        dag_api.execute_workflow()
    pass


run_workflow_dir(config_args['indir'])