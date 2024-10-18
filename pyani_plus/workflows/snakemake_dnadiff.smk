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
from pyani_plus.snakemake import snakemake_scheduler

indir_files = snakemake_scheduler.check_input_stems(config["indir"])


def get_genomeA(wildcards):
    return indir_files[wildcards.genomeA]


def get_genomeB(wildcards):
    return indir_files[wildcards.genomeB]


# The delta rule runs nucmer
# NOTE: We do not need to construct forward and reverse comparisons within this
# rule. This is done in the context of pyani_plus by specifying target files in
# the calling code
rule delta:
    params:
        nucmer=config["nucmer"],
        indir=config["indir"],
        outdir=config["outdir"],
    input:
        genomeA=get_genomeA,
        genomeB=get_genomeB,
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.delta",
    shell:
        "chronic {params.nucmer} -p {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB} --maxmatch {input.genomeB} {input.genomeA}"


# The filter rule runs delta-filter wrapper for nucmer
# NOTE: This rule is used for dnadiff in the context of pyani_plus. The .filter file
# is used to generate delta-filter files, which are used to replicate AlignedBases.
rule filter:
    params:
        delta_filter=config["delta_filter"],
    input:
        "{outdir}/{genomeA}_vs_{genomeB}.delta",
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.filter",
    shell:
        # Do not use chronic where we want stdout captured to file
        "{params.delta_filter} -m {input} > {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.filter"


# The rule runs both show-diff and show-coords wrappers for nucmer
# NOTE: This rule is used for dnadiff in the context of pyani_plus. The .filter file
# is used to generate both show-diff and show-coords files, which are used to
# replicate AlignedBases, identity etc.
rule show_diff_and_coords:
    params:
        db=config["db"],
        show_diff=config["show_diff"],
        show_coords=config["show_coords"],
        outdir=config["outdir"],
    input:
        "{outdir}/{genomeA}_vs_{genomeB}.filter",
        genomeA=get_genomeA,
        genomeB=get_genomeB,
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.qdiff",
        "{outdir}/{genomeA}_vs_{genomeB}.mcoords",
    shell:
        """
        {params.show_diff} -qH {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.filter > {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.qdiff
        {params.show_coords} -rclTH {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.filter > {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.mcoords &&
        .pyani-plus-private-cli log-dnadiff --database {params.db} \
            --query-fasta {input.genomeA} --subject-fasta {input.genomeB} \
            --mcoords {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.mcoords --qdiff {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.qdiff
        """
