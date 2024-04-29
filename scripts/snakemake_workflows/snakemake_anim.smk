
# The nucmer rule runs nucmer 
# NOTE: We do not need to construct forward and reverse comparisions separately
# This is done by secifying target files in `run_anim_workflow.py` and rule all.
rule nucmer:
    params:
        indir = config['indir'],
        mode = config['mode'],
        outdir = config['outdir']
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