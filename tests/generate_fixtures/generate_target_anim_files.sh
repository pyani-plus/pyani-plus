#!/bin/bash

# This bash script generates target files for anim comparisions.
# Genomes are compared in both directions (forward and reverse)
# using nucmer and delta-filter with --mum parameter (nucmer)
# and -1 parameter (delta-filter). 

# Paths to directories (eg, input sequences, delta and filter)
input_dir=../fixtures/sequences
delta_dir=../fixtures/anim/targets/delta
filter_dir=../fixtures/anim/targets/filter

#Running comparisions
nucmer -p $delta_dir/NC_002696_vs_NC_010338 --mum $input_dir/NC_002696.fna $input_dir/NC_010338.fna
delta-filter -1 $delta_dir/NC_002696_vs_NC_010338.delta > $filter_dir/NC_002696_vs_NC_010338.filter

nucmer -p $delta_dir/NC_002696_vs_NC_011916 --mum $input_dir/NC_002696.fna $input_dir/NC_011916.fna
delta-filter -1 $delta_dir/NC_002696_vs_NC_011916.delta > $filter_dir/NC_002696_vs_NC_011916.filter 

nucmer -p $delta_dir/NC_002696_vs_NC_014100 --mum $input_dir/NC_002696.fna $input_dir/NC_014100.fna
delta-filter -1 $delta_dir/NC_002696_vs_NC_014100.delta > $filter_dir/NC_002696_vs_NC_014100.filter 

nucmer -p $delta_dir/NC_010338_vs_NC_002696 --mum $input_dir/NC_010338.fna $input_dir/NC_002696.fna
delta-filter -1 $delta_dir/NC_010338_vs_NC_002696.delta > $filter_dir/NC_010338_vs_NC_002696.filter 

nucmer -p $delta_dir/NC_010338_vs_NC_011916 --mum $input_dir/NC_010338.fna $input_dir/NC_011916.fna
delta-filter -1 $delta_dir/NC_010338_vs_NC_011916.delta > $filter_dir/NC_010338_vs_NC_011916.filter 

nucmer -p $delta_dir/NC_010338_vs_NC_014100 --mum $input_dir/NC_010338.fna $input_dir/NC_014100.fna
delta-filter -1 $delta_dir/NC_010338_vs_NC_014100.delta > $filter_dir/NC_010338_vs_NC_014100.filter

nucmer -p $delta_dir/NC_011916_vs_NC_002696 --mum $input_dir/NC_011916.fna $input_dir/NC_002696.fna
delta-filter -1 $delta_dir/NC_011916_vs_NC_002696.delta > $filter_dir/NC_011916_vs_NC_002696.filter

nucmer -p $delta_dir/NC_011916_vs_NC_010338 --mum $input_dir/NC_011916.fna $input_dir/NC_010338.fna
delta-filter -1 $delta_dir/NC_011916_vs_NC_010338.delta > $filter_dir/NC_011916_vs_NC_010338.filter

nucmer -p $delta_dir/NC_011916_vs_NC_014100 --mum $input_dir/NC_011916.fna $input_dir/NC_014100.fna
delta-filter -1 $delta_dir/NC_011916_vs_NC_014100.delta > $filter_dir/NC_011916_vs_NC_014100.filter

nucmer -p $delta_dir/NC_014100_vs_NC_002696 --mum $input_dir/NC_014100.fna $input_dir/NC_002696.fna
delta-filter -1 $delta_dir/NC_014100_vs_NC_002696.delta > $filter_dir/NC_014100_vs_NC_002696.filter

nucmer -p $delta_dir/NC_014100_vs_NC_010338 --mum $input_dir/NC_014100.fna $input_dir/NC_010338.fna
delta-filter -1 $delta_dir/NC_014100_vs_NC_010338.delta > $filter_dir/NC_014100_vs_NC_010338.filter

nucmer -p $delta_dir/NC_014100_vs_NC_011916 --mum $input_dir/NC_014100.fna $input_dir/NC_011916.fna
delta-filter -1 $delta_dir/NC_014100_vs_NC_011916.delta > $filter_dir/NC_014100_vs_NC_011916.filter