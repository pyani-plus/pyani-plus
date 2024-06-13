# nucmer runs with the --maxmatch parameter to find initial matches
# regardless of their uniqueness.
# delta-filter runs with the -m parameter to filter only M-to-M
# matches.
# show-diff runs with the -q parameter to obtain structural
# differences for each sequence in the reference.
# NOTE: Depending on what values we wish to report using dnadiff
# subcommand, we might need to generate show-diff files for all
# sequences in the query too.


# Imports
import os

from pathlib import Path
from itertools import permutations

# Paths to directories (eg. input sequences, outputs for delta, filter...)
INPUT_DIR = "../fixtures/sequences"
DELTA_DIR = "../fixtures/dnadiff/targets/delta"
FILTER_DIR = "../fixtures/dnadiff/targets/filter"
SHOW_DIFF_DIR = "../fixtures/dnadiff/targets/show_diff"
SHOW_COORDS_DIR = "../fixtures/dnadiff/targets/show_coords"

# Running comparisions
comparisions = permutations([_.stem for _ in Path(INPUT_DIR).glob("*")], 2)
inputs = {_.stem: _ for _ in Path(INPUT_DIR).glob("*")}

for genomes in comparisions:
    os.system(
        f"nucmer -p {DELTA_DIR + '/' + '_vs_'.join(genomes)} --maxmatch {inputs[genomes[0]]} {inputs[genomes[1]]}"
    )
    os.system(
        f"delta-filter -m {DELTA_DIR + '/' + '_vs_'.join(genomes)}.delta > {FILTER_DIR + '/' + '_vs_'.join(genomes)}.filter"
    )
    os.system(
        f"show-diff -rH {FILTER_DIR + '/' + '_vs_'.join(genomes)}.filter > {SHOW_DIFF_DIR + '/' + '_vs_'.join(genomes)}.rdiff"
    )
    os.system(
        f"show-coords {FILTER_DIR + '/' + '_vs_'.join(genomes)}.filter > {SHOW_COORDS_DIR + '/' + '_vs_'.join(genomes)}.mcoords"
    )
