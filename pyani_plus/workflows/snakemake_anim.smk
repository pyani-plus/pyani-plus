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


# The ANIm rule runs nucmer and delta-filter
# NOTE: We do not need to construct forward and reverse comparisons within this
# rule. This is done in the context of pyani_plus by specifying target files in
# the calling code.
# Here mode will be "mum" (default) or "maxmatch", meaning nucmer --mum etc.
# The constant -1 option is used for 1-to-1 alignments in the delta-filter,
# with no other options available for the end user.
rule ANIm:
    params:
        db=config["db"],
        run_id=config["run_id"],
        mode=config["mode"],
        nucmer=config["nucmer"],
        delta_filter=config["delta_filter"],
        indir=config["indir"],
        outdir=config["outdir"],
    input:
        genomeA=get_genomeA,
        genomeB=get_genomeB,
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.filter",
    shell:
        """
        {params.nucmer} -p {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB} \
            --{params.mode} {input.genomeB} {input.genomeA} \
            > {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}_nucmer.log 2>&1 &&
        {params.delta_filter} \
            -1 {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.delta \
             > {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.filter &&
        .pyani-plus-private-cli log-anim --quiet \
            --database {params.db} --run-id {params.run_id} \
            --query-fasta {input.genomeA} --subject-fasta {input.genomeB} \
            --deltafilter {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.filter
        """
