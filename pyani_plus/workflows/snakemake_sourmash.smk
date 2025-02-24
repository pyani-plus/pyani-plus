# The MIT License
#
# Copyright (c) 2024-2025 University of Strathclyde
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""Snakemake workflow for sourmash"""


def get_fasta(wildcards):
    return str(config["md5_to_filename"][wildcards.genome_hash])


# The sketch rule runs the sourmash branchwater equivalent to "sourmash sketch"
rule sketch:
    params:
        outdir=config["outdir"],
        extra=config["extra"],  # This will consist of either `scaled=X` or `num=X`.
        kmersize=config["kmersize"],
    input:
        get_fasta,
    output:
        "{outdir}/{genome_hash}.sig",
    shell:
        """
        sourmash scripts singlesketch -I DNA \
            -p 'k={params.kmersize},{params.extra}' \
            -n "{wildcards.genome_hash}" \
            {input:q} -o "{output}" > "{output}.log" 2>&1
        """


# The compare rule runs the branchwater equivalent of "sourmash compare"
rule compare:
    params:
        db=config["db"],
        run_id=config["run_id"],
        outdir=config["outdir"],
        kmersize=config["kmersize"],
        temp=config["temp"],
        extra=config["extra"],  # This will consist of either `scaled=X` or `num=X`.
    input:
        expand("{{outdir}}/{genome}.sig", genome=sorted(config["md5_to_filename"])),
    output:
        "{outdir}/manysearch.csv",
    shell:
        """
        sourmash sig collect --quiet -F csv -o {wildcards.outdir}/all_sigs.csv {input:q} > {output}.log 2>&1 &&
        sourmash scripts manysearch -m DNA --quiet -t 0 -o {output} {wildcards.outdir}/all_sigs.csv {wildcards.outdir}/all_sigs.csv >> {output}.log 2>&1 &&
        .pyani-plus-private-cli log-sourmash --quiet \
            --database "{params.db}" --run-id {params.run_id} \
            --manysearch "{output}"
        """
