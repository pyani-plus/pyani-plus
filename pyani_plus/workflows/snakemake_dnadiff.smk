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
rule delta:
    params:
        indir=config["indir"],
        outdir=config["outdir"],
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.delta",
    input:
        genomeA=get_genomeA,
        genomeB=get_genomeB
    run:
        shell(
            "nucmer -p {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB} --maxmatch {input.genomeA} {input.genomeB}"
        )


# The filter rule runs delta-filter wrapper for nucmer
# NOTE: This rule is used for dnadiff in the context of pyani_plus. The .filter file
# is used to generate delta-filter files, which are used to replicate AlignedBases. 
rule filter:
    input:
        "{outdir}/{genomeA}_vs_{genomeB}.delta",
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.filter",
    run:
        shell(
            "delta-filter -m {input} > {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.filter"
        )

# The filter rule runs show-diff wrapper for nucmer
# NOTE: This rule is used for dnadiff in the context of pyani_plus. The .filter file
# is used to generate show-diff files, which are used to replicate AlignedBases. 
rule show_diff:
    input:
        "{outdir}/{genomeA}_vs_{genomeB}.filter",
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.rdiff",
    run:
        shell(
            "show-diff -rH {input} > {output}"
        )

# The filter rule runs show-coords wrapper for nucmer
# NOTE: This rule is used for dnadiff in the context of pyani_plus. The .filter file
# is used to generate show-coords files, which are used to replicate AlignedBases. 
rule show_coords:
    input:
        "{outdir}/{genomeA}_vs_{genomeB}.filter",
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.mcoords",
    run:
        shell(
            "show-coords {input} > {output}"
        )
