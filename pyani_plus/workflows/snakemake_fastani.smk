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
