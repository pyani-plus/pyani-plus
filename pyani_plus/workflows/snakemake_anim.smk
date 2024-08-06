from pyani_plus.snakemake import snakemake_scheduler

indir_files = snakemake_scheduler.check_input_stems(config['indir'])

def get_genomeA(wildcards):
    return indir_files[wildcards.genomeA]

def get_genomeB(wildcards):
    return indir_files[wildcards.genomeB]




# The delta rule runs nucmer
# NOTE: We do not need to construct forward and reverse comparisions within this
# rule. This is done in the context of pyani_plus by specifying target files in
# the calling code

from pyani_plus.snakemake import anim

indir_files = anim.check_input_stems(config['indir'])

def get_genomeA(wildcards):
    return indir_files[wildcards.genomeA]

def get_genomeB(wildcards):
    return indir_files[wildcards.genomeB]


rule delta:
    params:
        indir=config["indir"],
        mode=config["mode"],
        outdir=config["outdir"]
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.delta",
    input:
        genomeA=get_genomeA,
        genomeB=get_genomeB
    run:
        shell(
            "nucmer -p {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB} --{params.mode} {input.genomeA} {input.genomeB}"
        )


# The filter rule runs delta-filter wrapper for nucmer
# NOTE: This rule is used for ANIm in the context of pyani_plus. The .filter file
# is used to calculate ANI values
rule filter:
    input:
        "{outdir}/{genomeA}_vs_{genomeB}.delta",
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.filter",
    run:
        shell(
            "delta-filter -1 {input} > {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.filter"
        )
