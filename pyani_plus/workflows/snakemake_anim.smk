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
"""Snakemake workflow for ANIm"""
from pyani_plus.workflows import check_input_stems

indir_files = check_input_stems(config["indir"])


def get_genomeA(wildcards):
    return indir_files[wildcards.genomeA]


def get_genomeB(wildcards):
    return indir_files[wildcards.genomeB]


# The delta rule runs nucmer
# NOTE: We do not need to construct forward and reverse comparisons within this
# rule. This is done in the context of pyani_plus by specifying target files in
# the calling code
# Here mode will be "mum" (default) or "maxmatch", meaning nucmer --mum etc.
rule delta:
    params:
        nucmer=config["nucmer"],
        mode=config["mode"],
        indir=config["indir"],
        outdir=config["outdir"],
    input:
        genomeA=get_genomeA,
        genomeB=get_genomeB,
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.delta",
    shell:
        """
        {params.nucmer} -p {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB} \
            --{params.mode} {input.genomeB} {input.genomeA} > {output}.log 2>&1
        """


# The filter rule runs delta-filter wrapper for nucmer
# NOTE: This rule is used for ANIm in the context of pyani_plus. The .filter file
# is used to calculate ANI values
rule filter:
    params:
        db=config["db"],
        mode=config["mode"],
        delta_filter=config["delta_filter"],
        outdir=config["outdir"],
    input:
        "{outdir}/{genomeA}_vs_{genomeB}.delta",
        genomeA=get_genomeA,
        genomeB=get_genomeB,
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.filter",
    shell:
        """
        {params.delta_filter} \
            -1 {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.delta \
             > {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.filter &&
        .pyani-plus-private-cli log-anim --quiet --database {params.db} \
            --query-fasta {input.genomeA} --subject-fasta {input.genomeB} \
            --deltafilter {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.filter \
            --mode {params.mode}
        """
