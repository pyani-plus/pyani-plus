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

# This returns a dict mapping stems to full paths
indir_files = check_input_stems(config["indir"])


def get_genomeB(wildcards):
    return indir_files[wildcards.genomeB]


# Need a file with a single FASTA file per line, one for each query.
# This could be all the FASTA files in the input directory, but can
# be a subset - consider use case expanding a DB from N genomes to
# N+1 genomes.
rule genomes_list:
    params:
        db=config["db"],
        run_id=config["run_id"],
        indir=config["indir"],
    input:
        genomeB=get_genomeB,
    output:
        "{outdir}/genome_list_for_{genomeB}.txt",
    shell:
        """
        .pyani-plus-private-cli build-query-list --quiet --self \
            --database {params.db} --run-id {params.run_id} \
            --subject {input} --fasta {params.indir} > {output}
        """


# The fastani rule runs fastANI doing all queries vs one subject (reference)
rule fastani:
    params:
        db=config["db"],
        run_id=config["run_id"],
        fastani=config["fastani"],
        indir=config["indir"],
        outdir=config["outdir"],
        fragsize=config["fragsize"],
        kmersize=config["kmersize"],
        minmatch=config["minmatch"],
    input:
        queries="{outdir}/genome_list_for_{genomeB}.txt",
        genomeB=get_genomeB,
    output:
        "{outdir}/all_vs_{genomeB}.fastani",
    shell:
        """
        {params.fastani} --ql "{input.queries}" -r {input.genomeB} \
            -o {output} --fragLen {params.fragsize} -k {params.kmersize} \
            --minFraction {params.minmatch} > {output}.log 2>&1 &&
        .pyani-plus-private-cli log-fastani --quiet \
            --database {params.db} --run-id {params.run_id} \
            --fastani {output}
        """
