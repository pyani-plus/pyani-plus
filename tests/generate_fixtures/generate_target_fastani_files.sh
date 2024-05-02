#!/bin/env bash

# This bash script generates target files for fastani comparisions.
# Genomes are compared in both directions (forward and reverse)
# using fastANI.

# Additional parameters (eg, input sequences, fastANI outputs, k-mer sizes...)
INPUT_DIR=../fixtures/sequences
FASTANI_DIR=../fixtures/fastani/targets
FRAG_LEN=3000
KMER_SIZE=16
MIN_FRAC=0.2


#Running comparisions
fastANI -q $INPUT_DIR/NC_002696.fna -r $INPUT_DIR/NC_010338.fna -o $FASTANI_DIR/NC_002696_vs_NC_010338.fastani --fragLen $FRAG_LEN -k $KMER_SIZE --minFraction $MIN_FRAC

fastANI -q $INPUT_DIR/NC_002696.fna -r $INPUT_DIR/NC_011916.fna -o $FASTANI_DIR/NC_002696_vs_NC_011916.fastani --fragLen $FRAG_LEN -k $KMER_SIZE --minFraction $MIN_FRAC

fastANI -q $INPUT_DIR/NC_002696.fna -r $INPUT_DIR/NC_014100.fna -o $FASTANI_DIR/NC_002696_vs_NC_014100.fastani --fragLen $FRAG_LEN -k $KMER_SIZE --minFraction $MIN_FRAC

fastANI -q $INPUT_DIR/NC_010338.fna -r $INPUT_DIR/NC_002696.fna -o $FASTANI_DIR/NC_010338_vs_NC_002696.fastani --fragLen $FRAG_LEN -k $KMER_SIZE --minFraction $MIN_FRAC

fastANI -q $INPUT_DIR/NC_010338.fna -r $INPUT_DIR/NC_011916.fna -o $FASTANI_DIR/NC_010338_vs_NC_011916.fastani --fragLen $FRAG_LEN -k $KMER_SIZE --minFraction $MIN_FRAC

fastANI -q $INPUT_DIR/NC_010338.fna -r $INPUT_DIR/NC_014100.fna -o $FASTANI_DIR/NC_010338_vs_NC_014100.fastani --fragLen $FRAG_LEN -k $KMER_SIZE --minFraction $MIN_FRAC

fastANI -q $INPUT_DIR/NC_011916.fna -r $INPUT_DIR/NC_002696.fna -o $FASTANI_DIR/NC_011916_vs_NC_002696.fastani --fragLen $FRAG_LEN -k $KMER_SIZE --minFraction $MIN_FRAC

fastANI -q $INPUT_DIR/NC_011916.fna -r $INPUT_DIR/NC_010338.fna -o $FASTANI_DIR/NC_011916_vs_NC_010338.fastani --fragLen $FRAG_LEN -k $KMER_SIZE --minFraction $MIN_FRAC

fastANI -q $INPUT_DIR/NC_011916.fna -r $INPUT_DIR/NC_014100.fna -o $FASTANI_DIR/NC_011916_vs_NC_014100.fastani --fragLen $FRAG_LEN -k $KMER_SIZE --minFraction $MIN_FRAC

fastANI -q $INPUT_DIR/NC_014100.fna -r $INPUT_DIR/NC_002696.fna -o $FASTANI_DIR/NC_014100_vs_NC_002696.fastani --fragLen $FRAG_LEN -k $KMER_SIZE --minFraction $MIN_FRAC

fastANI -q $INPUT_DIR/NC_014100.fna -r $INPUT_DIR/NC_010338.fna -o $FASTANI_DIR/NC_014100_vs_NC_010338.fastani --fragLen $FRAG_LEN -k $KMER_SIZE --minFraction $MIN_FRAC

fastANI -q $INPUT_DIR/NC_014100.fna -r $INPUT_DIR/NC_011916.fna -o $FASTANI_DIR/NC_014100_vs_NC_011916.fastani --fragLen $FRAG_LEN -k $KMER_SIZE --minFraction $MIN_FRAC
