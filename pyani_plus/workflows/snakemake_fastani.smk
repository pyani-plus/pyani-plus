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
from pyani_plus.workflows import check_input_stems

indir_files = check_input_stems(config["indir"])


def get_genomeA(wildcards):
    return indir_files[wildcards.genomeA]


def get_genomeB(wildcards):
    return indir_files[wildcards.genomeB]


# The fastani rule runs fastANI
rule fastani:
    params:
        db=config["db"],
        fastani=config["fastani"],
        indir=config["indir"],
        fragsize=config["fragsize"],
        kmersize=config["kmersize"],
        minmatch=config["minmatch"],
    input:
        genomeA=get_genomeA,
        genomeB=get_genomeB,
    output:
        "{outdir}/{genomeA}_vs_{genomeB}.fastani",
    shell:
        """
        chronic {params.fastani} -q {input.genomeA} -r {input.genomeB} \
            -o {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.fastani \
            --fragLen {params.fragsize} -k {params.kmersize} --minFraction {params.minmatch} &&
        chronic .pyani-plus-private-cli log-fastani --database {params.db} \
            --query-fasta {input.genomeA} --subject-fasta {input.genomeB} \
            --fastani {wildcards.outdir}/{wildcards.genomeA}_vs_{wildcards.genomeB}.fastani \
            --fragsize {params.fragsize} --kmersize {params.kmersize} --minmatch {params.minmatch}
        """
