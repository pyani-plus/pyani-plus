from pyani_plus.snakemake import snakemake_scheduler

indir_files = snakemake_scheduler.check_input_stems(config["indir"])


def get_genomeA(wildcards):
    return indir_files[wildcards.genomeA]


def get_genomeB(wildcards):
    return indir_files[wildcards.genomeB]


# The fastani rule runs fastANI
rule fastani:
    params:
        fastani=config["fastani"],
        indir=config["indir"],
        fragLen=config["fragLen"],
        kmerSize=config["kmerSize"],
        minFrac=config["minFrac"],
    input:
        genomeA=get_genomeA,
        genomeB=get_genomeB,
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.fastani",
    shell:
        "{params.fastani} -q {input.genomeA} -r {input.genomeB} -o {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.fastani --fragLen {params.fragLen} -k {params.kmerSize} --minFraction {params.minFrac}"
