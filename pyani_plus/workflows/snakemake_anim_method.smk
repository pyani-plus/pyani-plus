rule anim_method:
    params:
        indir=config["indir"],
        outdir=config["outdir"],
    output:
        "{outdir}/{comp}.txt",
    run:
        shell(
            "python -m pyani_plus.method.method_anim parse_delta {params.indir}/{wildcards.comp}.filter {output}"
        )
