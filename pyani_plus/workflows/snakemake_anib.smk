# The MIT License
#
# Copyright (c) 2024 University of Strathclyde
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
"""Snakemake workflow for ANIb"""
from pyani_plus.workflows import check_input_stems

indir_files = check_input_stems(config["indir"])


def get_genomeA(wildcards):
    return indir_files[wildcards.genomeA]


def get_genomeB(wildcards):
    return indir_files[wildcards.genomeB]


# The fragments FASTA file is used for the ANIb query.
# Should we put use fragment size in the filename?
rule fragment:
    params:
        fragsize=config["fragsize"],
        indir=config["indir"],
        outdir=config["outdir"],
    input:
        genomeA=get_genomeA,
    output:
        "{outdir}/{genomeA}-fragments.fna",
    shell:
        ".pyani-plus-private-cli fragment-fasta --quiet --outdir {wildcards.outdir} --fragsize {params.fragsize} {input.genomeA}"


# The nucleotide database of the FASTA file is used for the ANIb reference.
# A nucleotide BLAST database is a set of files <stem>.n*
# Picking <stem>.njs as the one we'll track for snakemake
# as this is plain text, specifically a JSON summary
rule blastdb:
    params:
        makeblastdb=config["makeblastdb"],
        indir=config["indir"],
        outdir=config["outdir"],
    input:
        genomeB=get_genomeB,
    output:
        "{outdir}/{genomeB}.njs",
    shell:
        """
        {params.makeblastdb} -in {input.genomeB} -input_type fasta -dbtype nucl \
            -title {wildcards.genomeB} -out {wildcards.outdir}/{wildcards.genomeB} \
            > {wildcards.outdir}/{wildcards.genomeB}.log
        """


# For ANIb query the fragments FASTA from genomeA against a BLAST DB of genomeB.
rule blastn:
    params:
        db=config["db"],
        run_id=config["run_id"],
        blastn=config["blastn"],
        fragsize=config["fragsize"],
        indir=config["indir"],
        outdir=config["outdir"],
    input:
        "{outdir}/{genomeA}-fragments.fna",
        "{outdir}/{genomeB}.njs",
        genomeA=get_genomeA,
        genomeB=get_genomeB,
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.tsv",
    shell:
        """
        {params.blastn} -query {wildcards.outdir}/{wildcards.genomeA}-fragments.fna \
            -db {wildcards.outdir}/{wildcards.genomeB} -out {output} -task blastn \
            -outfmt '6 qseqid sseqid length mismatch pident nident qlen slen \
                     qstart qend sstart send positive ppos gaps' \
            -xdrop_gap_final 150 -dust no -evalue 1e-15 > {output}.log &&
        .pyani-plus-private-cli log-anib --quiet \
            --database {params.db} --run-id {params.run_id} \
            --query-fasta {input.genomeA} --subject-fasta {input.genomeB} \
            --blastn {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.tsv
        """
