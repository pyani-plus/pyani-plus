#!/bin/env bash
#
# This bash script generates target files for anim comparisions.
# Genomes are compared in both directions (forward and reverse)
# using nucmer and delta-filter.
#
# nucmer runs with the --mum parameter to find initial maximal
# unique matches.
# delta-filter runs with the -1 parameter to filter only 1-to-1
# matches

# Paths to directories (eg, input sequences, delta and filter)
INPUT_DIR=../fixtures/sequences
DELTA_DIR=../fixtures/anim/targets/delta
FILTER_DIR=../fixtures/anim/targets/filter

#Running comparisions
nucmer -p $DELTA_DIR/NC_002696_vs_NC_010338 --mum $INPUT_DIR/NC_002696.fna $INPUT_DIR/NC_010338.fna
delta-filter -1 $DELTA_DIR/NC_002696_vs_NC_010338.delta > $FILTER_DIR/NC_002696_vs_NC_010338.filter

nucmer -p $DELTA_DIR/NC_002696_vs_NC_011916 --mum $INPUT_DIR/NC_002696.fna $INPUT_DIR/NC_011916.fna
delta-filter -1 $DELTA_DIR/NC_002696_vs_NC_011916.delta > $FILTER_DIR/NC_002696_vs_NC_011916.filter 

nucmer -p $DELTA_DIR/NC_002696_vs_NC_014100 --mum $INPUT_DIR/NC_002696.fna $INPUT_DIR/NC_014100.fna
delta-filter -1 $DELTA_DIR/NC_002696_vs_NC_014100.delta > $FILTER_DIR/NC_002696_vs_NC_014100.filter 

nucmer -p $DELTA_DIR/NC_010338_vs_NC_002696 --mum $INPUT_DIR/NC_010338.fna $INPUT_DIR/NC_002696.fna
delta-filter -1 $DELTA_DIR/NC_010338_vs_NC_002696.delta > $FILTER_DIR/NC_010338_vs_NC_002696.filter 

nucmer -p $DELTA_DIR/NC_010338_vs_NC_011916 --mum $INPUT_DIR/NC_010338.fna $INPUT_DIR/NC_011916.fna
delta-filter -1 $DELTA_DIR/NC_010338_vs_NC_011916.delta > $FILTER_DIR/NC_010338_vs_NC_011916.filter 

nucmer -p $DELTA_DIR/NC_010338_vs_NC_014100 --mum $INPUT_DIR/NC_010338.fna $INPUT_DIR/NC_014100.fna
delta-filter -1 $DELTA_DIR/NC_010338_vs_NC_014100.delta > $FILTER_DIR/NC_010338_vs_NC_014100.filter

nucmer -p $DELTA_DIR/NC_011916_vs_NC_002696 --mum $INPUT_DIR/NC_011916.fna $INPUT_DIR/NC_002696.fna
delta-filter -1 $DELTA_DIR/NC_011916_vs_NC_002696.delta > $FILTER_DIR/NC_011916_vs_NC_002696.filter

nucmer -p $DELTA_DIR/NC_011916_vs_NC_010338 --mum $INPUT_DIR/NC_011916.fna $INPUT_DIR/NC_010338.fna
delta-filter -1 $DELTA_DIR/NC_011916_vs_NC_010338.delta > $FILTER_DIR/NC_011916_vs_NC_010338.filter

nucmer -p $DELTA_DIR/NC_011916_vs_NC_014100 --mum $INPUT_DIR/NC_011916.fna $INPUT_DIR/NC_014100.fna
delta-filter -1 $DELTA_DIR/NC_011916_vs_NC_014100.delta > $FILTER_DIR/NC_011916_vs_NC_014100.filter

nucmer -p $DELTA_DIR/NC_014100_vs_NC_002696 --mum $INPUT_DIR/NC_014100.fna $INPUT_DIR/NC_002696.fna
delta-filter -1 $DELTA_DIR/NC_014100_vs_NC_002696.delta > $FILTER_DIR/NC_014100_vs_NC_002696.filter

nucmer -p $DELTA_DIR/NC_014100_vs_NC_010338 --mum $INPUT_DIR/NC_014100.fna $INPUT_DIR/NC_010338.fna
delta-filter -1 $DELTA_DIR/NC_014100_vs_NC_010338.delta > $FILTER_DIR/NC_014100_vs_NC_010338.filter

nucmer -p $DELTA_DIR/NC_014100_vs_NC_011916 --mum $INPUT_DIR/NC_014100.fna $INPUT_DIR/NC_011916.fna
delta-filter -1 $DELTA_DIR/NC_014100_vs_NC_011916.delta > $FILTER_DIR/NC_014100_vs_NC_011916.filter