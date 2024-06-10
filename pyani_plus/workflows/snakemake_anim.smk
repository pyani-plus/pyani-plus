# The delta rule runs nucmer
# NOTE: We do not need to construct forward and reverse comparisions within this
# rule. This is done in the context of pyani_plus by specifying target files in
# the calling code

# Define acceptable extensions
extensions = {".fna", ".fasta", ".fas"}

def get_correct_file(genome):
    "Return file for a given genome with acceptable extention."

    for ext in extensions:
        file = Path(config["indir"]) / f"{genome}{ext}"
        if file.exists():
            return file

def resolve_genomeA(wildcards):
    "Return genome A."
    return get_correct_file(wildcards.genomeA)

def resolve_genomeB(wildcards):
    "Return genome B."
    return get_correct_file(wildcards.genomeB)


rule delta:
    params:
        indir=config["indir"],
        mode=config["mode"],
        outdir=config["outdir"]
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.delta"
    input:
        genomeA=resolve_genomeA,
        genomeB=resolve_genomeB
    shell:
        "nucmer -p {params.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB} --{params.mode} {input.genomeA} {input.genomeB}"


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
