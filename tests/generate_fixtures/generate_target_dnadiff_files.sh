#!/bin/env bash
#
# This bash script generates target files for anim comparisions.
# Genomes are compared in both directions (forward and reverse)
# using nucmer and delta-filter.
#
# nucmer runs with the --maxmatch parameter to find initial matches
# regardless of their uniqueness.
# delta-filter runs with the -m parameter to filter only M-to-M
# matches
# show-diff runs with the -q parameter to obtain structural 
# differences for each sequence in the reference. 
# NOTE: Depending on what values we wish to report using dnadiff 
# subcommand, we might need to generate show-diff files for all 
# sequences in the query in the future. 



# Paths to directories (eg, input sequences, delta and filter)
INPUT_DIR=../fixtures/sequences
DELTA_DIR=../fixtures/dnadiff/targets/delta
FILTER_DIR=../fixtures/dnadiff/targets/filter
SHOW_DIFF_DIR=../fixtures/dnadiff/targets/show_diff

#Running comparisions
nucmer -p $DELTA_DIR/NC_002696_vs_NC_010338 --maxmatch $INPUT_DIR/NC_002696.fna $INPUT_DIR/NC_010338.fna
delta-filter -m $DELTA_DIR/NC_002696_vs_NC_010338.delta > $FILTER_DIR/NC_002696_vs_NC_010338.filter
show-diff -rH $FILTER_DIR/NC_002696_vs_NC_010338.filter > $SHOW_DIFF_DIR/NC_002696_vs_NC_010338.rdiff


nucmer -p $DELTA_DIR/NC_002696_vs_NC_011916 --maxmatch $INPUT_DIR/NC_002696.fna $INPUT_DIR/NC_011916.fna
delta-filter -m $DELTA_DIR/NC_002696_vs_NC_011916.delta > $FILTER_DIR/NC_002696_vs_NC_011916.filter
show-diff -rH $FILTER_DIR/NC_002696_vs_NC_011916.filter > $SHOW_DIFF_DIR/NC_002696_vs_NC_011916.rdiff


nucmer -p $DELTA_DIR/NC_002696_vs_NC_014100 --maxmatch $INPUT_DIR/NC_002696.fna $INPUT_DIR/NC_014100.fna
delta-filter -m $DELTA_DIR/NC_002696_vs_NC_014100.delta > $FILTER_DIR/NC_002696_vs_NC_014100.filter
show-diff -rH  $FILTER_DIR/NC_002696_vs_NC_014100.filter > $SHOW_DIFF_DIR/NC_002696_vs_NC_014100.rdiff


nucmer -p $DELTA_DIR/NC_010338_vs_NC_002696 --maxmatch $INPUT_DIR/NC_010338.fna $INPUT_DIR/NC_002696.fna
delta-filter -m $DELTA_DIR/NC_010338_vs_NC_002696.delta > $FILTER_DIR/NC_010338_vs_NC_002696.filter 
show-diff -rH  $FILTER_DIR/NC_010338_vs_NC_002696.filter > $SHOW_DIFF_DIR/NC_010338_vs_NC_002696.rdiff


nucmer -p $DELTA_DIR/NC_010338_vs_NC_011916 --maxmatch $INPUT_DIR/NC_010338.fna $INPUT_DIR/NC_011916.fna
delta-filter -m $DELTA_DIR/NC_010338_vs_NC_011916.delta > $FILTER_DIR/NC_010338_vs_NC_011916.filter 
show-diff -rH  $FILTER_DIR/NC_010338_vs_NC_011916.filter > $SHOW_DIFF_DIR/NC_010338_vs_NC_011916.rdiff


nucmer -p $DELTA_DIR/NC_010338_vs_NC_014100 --maxmatch $INPUT_DIR/NC_010338.fna $INPUT_DIR/NC_014100.fna
delta-filter -m $DELTA_DIR/NC_010338_vs_NC_014100.delta > $FILTER_DIR/NC_010338_vs_NC_014100.filter
show-diff -rH  $FILTER_DIR/NC_010338_vs_NC_014100.filter > $SHOW_DIFF_DIR/NC_010338_vs_NC_014100.rdiff


nucmer -p $DELTA_DIR/NC_011916_vs_NC_002696 --maxmatch $INPUT_DIR/NC_011916.fna $INPUT_DIR/NC_002696.fna
delta-filter -m $DELTA_DIR/NC_011916_vs_NC_002696.delta > $FILTER_DIR/NC_011916_vs_NC_002696.filter
show-diff -rH  $FILTER_DIR/NC_011916_vs_NC_002696.filter > $SHOW_DIFF_DIR/NC_011916_vs_NC_002696.rdiff


nucmer -p $DELTA_DIR/NC_011916_vs_NC_010338 --maxmatch $INPUT_DIR/NC_011916.fna $INPUT_DIR/NC_010338.fna
delta-filter -m $DELTA_DIR/NC_011916_vs_NC_010338.delta > $FILTER_DIR/NC_011916_vs_NC_010338.filter
show-diff -rH  $FILTER_DIR/NC_011916_vs_NC_010338.filter > $SHOW_DIFF_DIR/NC_011916_vs_NC_010338.rdiff


nucmer -p $DELTA_DIR/NC_011916_vs_NC_014100 --maxmatch $INPUT_DIR/NC_011916.fna $INPUT_DIR/NC_014100.fna
delta-filter -m $DELTA_DIR/NC_011916_vs_NC_014100.delta > $FILTER_DIR/NC_011916_vs_NC_014100.filter
show-diff -rH  $FILTER_DIR/NC_011916_vs_NC_014100.filter > $SHOW_DIFF_DIR/NC_011916_vs_NC_014100.rdiff


nucmer -p $DELTA_DIR/NC_014100_vs_NC_002696 --maxmatch $INPUT_DIR/NC_014100.fna $INPUT_DIR/NC_002696.fna
delta-filter -m $DELTA_DIR/NC_014100_vs_NC_002696.delta > $FILTER_DIR/NC_014100_vs_NC_002696.filter
show-diff -rH  $FILTER_DIR/NC_014100_vs_NC_002696.filter > $SHOW_DIFF_DIR/NC_014100_vs_NC_002696.rdiff


nucmer -p $DELTA_DIR/NC_014100_vs_NC_010338 --maxmatch $INPUT_DIR/NC_014100.fna $INPUT_DIR/NC_010338.fna
delta-filter -m $DELTA_DIR/NC_014100_vs_NC_010338.delta > $FILTER_DIR/NC_014100_vs_NC_010338.filter
show-diff -rH  $FILTER_DIR/NC_014100_vs_NC_010338.filter > $SHOW_DIFF_DIR/NC_014100_vs_NC_010338.rdiff


nucmer -p $DELTA_DIR/NC_014100_vs_NC_011916 --maxmatch $INPUT_DIR/NC_014100.fna $INPUT_DIR/NC_011916.fna
delta-filter -m $DELTA_DIR/NC_014100_vs_NC_011916.delta > $FILTER_DIR/NC_014100_vs_NC_011916.filter
show-diff -rH  $FILTER_DIR/NC_014100_vs_NC_011916.filter > $SHOW_DIFF_DIR/NC_014100_vs_NC_011916.rdiff