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
        "{params.nucmer} -p {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB} --maxmatch {input.genomeA} {input.genomeB}"


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
        "{params.delta_filter} -m {input} > {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.filter"


# The filter rule runs show-diff wrapper for nucmer
# NOTE: This rule is used for dnadiff in the context of pyani_plus. The .filter file
# is used to generate show-diff files, which are used to replicate AlignedBases.
rule show_diff:
    params:
        show_diff=config["show_diff"],
    input:
        "{outdir}/{genomeA}_vs_{genomeB}.filter",
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.qdiff",
    shell:
        "{params.show_diff} -qH {input} > {output}"


# The filter rule runs show-coords wrapper for nucmer
# NOTE: This rule is used for dnadiff in the context of pyani_plus. The .filter file
# is used to generate show-coords files, which are used to replicate AlignedBases.
rule show_coords:
    params:
        show_coords=config["show_coords"],
    input:
        "{outdir}/{genomeA}_vs_{genomeB}.filter",
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.mcoords",
    shell:
        "{params.show_coords} -rclTH {input} > {output}"
