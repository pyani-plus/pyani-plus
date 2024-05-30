# The fastani rule runs fastANI

rule fastani:
    params:
        indir=config["indir"],
        fragLen=config["fragLen"],
        kmerSize=config["kmerSize"],
        minFrac=config["minFrac"],


    output:
        "{outdir}/{genomeA}_vs_{genomeB}.fastani",
    run:
        shell(
            "fastANI -q {params.indir}/{wildcards.genomeA}.* -r {params.indir}/{wildcards.genomeB}.* -o {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.fastani --fragLen {params.fragLen} -k {params.kmerSize} --minFraction {params.minFrac}"
        )