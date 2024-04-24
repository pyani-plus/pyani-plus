# Example snakemake file performing pairwise comparisons
from itertools import permutations


# Pick up genome filestems
(GENOMES,) = glob_wildcards(str(config['indir'])+"/{genome}.fna")
CMPS = list(permutations(GENOMES, 2))  # all pairwise comparisons fwd and reverse


# Rule `all` defines all A vs B comparisons, the `nucmer` rule runs a
# single pairwise comparison at a time
# The `zip` argument to `expand()` prevents this function generating the
# product of every member of each list. Instead we have extracted each
# participant in all pairwise comparisons into separate lists
rule all:
    input:
        expand(
            "{outdir}/{genomeA}_vs_{genomeB}.filter",
            zip,
            genomeA=[_[0] for _ in CMPS],
            genomeB=[_[1] for _ in CMPS],
            outdir=config['outdir']
        ),


# The nucmer rule runs nucmer 
# NOTE: We do not need to construct forward and reverse comparisions separately
# This is done by secifying target files in `run_anim_workflow.py` and rule all.
rule nucmer:
    params:
        indir = config['indir'],
        mode = config['mode']
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.delta",
    run:
        shell(
            "nucmer -p {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB} --{params.mode} {params.indir}/{wildcards.genomeA}.fna {params.indir}/{wildcards.genomeB}.fna"
        )


#The filter rule runs delta-filter wraper for nucmer
# NOTE: by specifying delta files as input will enfore the order in wich the rules will be executed
# eg. nucmer rule to provide delta files, then filter rule to obtain filter files

rule filter:
    input:
        "{outdir}/{genomeA}_vs_{genomeB}.delta"
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.filter",
    run:
        shell(
            "delta-filter -1 {input} > {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.filter"
        )