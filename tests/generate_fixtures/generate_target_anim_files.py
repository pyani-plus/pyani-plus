# This bash script generates target files for anim comparisions.
# Genomes are compared in both directions (forward and reverse)
# using nucmer and delta-filter.
#
# nucmer runs with the --mum parameter to find initial maximal
# unique matches.
# delta-filter runs with the -1 parameter to filter only 1-to-1
# matches

# Imports
import os

from pathlib import Path
from itertools import permutations

# Paths to directories (eg, input sequences, delta and filter)
INPUT_DIR = "../fixtures/sequences"
DELTA_DIR = "../fixtures/anim/targets/delta"
FILTER_DIR = "../fixtures/anim/targets/filter"

# Running ANIm comparisions
comparisions = permutations([_.stem for _ in Path(INPUT_DIR).glob("*")], 2)
inputs = {_.stem: _ for _ in Path(INPUT_DIR).glob("*")}

for genomes in comparisions:
    os.system(
        f"nucmer -p {DELTA_DIR + '/' + '_vs_'.join(genomes)} --mum {inputs[genomes[0]]} {inputs[genomes[1]]}"
    )
    os.system(
        f"delta-filter -1 {DELTA_DIR + '/' + '_vs_'.join(genomes)}.delta > {FILTER_DIR + '/' + '_vs_'.join(genomes)}.filter"
    )
