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
"""Snakemake workflow for pyANI-plus compute methods split by column.

This assumes there are no job dependencies, or they have been resolved already.
"""
from pyani_plus.workflows import check_input_stems

indir_files = check_input_stems(config["indir"])


def get_genomeB(wildcards):
    return indir_files[wildcards.genomeB]


rule compute_column:
    params:
        db=config["db"],
        run_id=config["run_id"],
        outdir=config["outdir"],
        temp=config["temp"],
    input:
        genomeB=get_genomeB,
    output:
        "{outdir}/all_vs_{genomeB}.{method}",
    shell:
        """
        .pyani-plus-private-cli compute-column --quiet \
            --database "{params.db}" --run-id {params.run_id} \
            --subject "{input}" {params.temp} && touch "{output}"
        """
