# This bash script generates target files for fastani comparisions.
# Genomes are compared in both directions (forward and reverse)
# using fastANI.

# Imports
import os

from pathlib import Path
from itertools import permutations


# Parameters (eg, input sequences, fastANI outputs, k-mer sizes...)
INPUT_DIR = "../fixtures/sequences"
FASTANI_DIR = "../fixtures/fastani/targets"
FRAG_LEN = 3000
KMER_SIZE = 16
MIN_FRAC = 0.2

# Running comparisions
comparisions = permutations([_.stem for _ in Path(INPUT_DIR).glob("*")], 2)
inputs = {_.stem: _ for _ in Path(INPUT_DIR).glob("*")}

for genomes in comparisions:
    os.system(
        f"fastANI -q {inputs[genomes[0]]} -r {inputs[genomes[1]]} -o {FASTANI_DIR + '/' + '_vs_'.join(genomes)}.fastani --fragLen {FRAG_LEN} -k {KMER_SIZE} --minFraction {MIN_FRAC}"
    )
