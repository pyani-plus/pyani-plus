#!/bin/bash

# This bash script generates target files for fastani comparisions.
# Genomes are compared in both directions (forward and reverse)
# using fastANI.

# Paths to directories (eg, input sequences, fastANI outputs)
input_dir=../fixtures/sequences
fastani_dir=../fixtures/fastani/targets
fragLen=3000
kmerSize=16
minFrac=0.2


#Running comparisions
fastANI -q $input_dir/NC_002696.fna -r $input_dir/NC_010338.fna -o $fastani_dir/NC_002696_vs_NC_010338.fastani --fragLen $fragLen -k $kmerSize --minFraction $minFrac


fastANI -q $input_dir/NC_002696.fna -r $input_dir/NC_011916.fna -o $fastani_dir/NC_002696_vs_NC_011916.fastani --fragLen $fragLen -k $kmerSize --minFraction $minFrac


fastANI -q $input_dir/NC_002696.fna -r $input_dir/NC_014100.fna -o $fastani_dir/NC_002696_vs_NC_014100.fastani --fragLen $fragLen -k $kmerSize --minFraction $minFrac


fastANI -q $input_dir/NC_010338.fna -r $input_dir/NC_002696.fna -o $fastani_dir/NC_010338_vs_NC_002696.fastani --fragLen $fragLen -k $kmerSize --minFraction $minFrac


fastANI -q $input_dir/NC_010338.fna -r $input_dir/NC_011916.fna -o $fastani_dir/NC_010338_vs_NC_011916.fastani --fragLen $fragLen -k $kmerSize --minFraction $minFrac


fastANI -q $input_dir/NC_010338.fna -r $input_dir/NC_014100.fna -o $fastani_dir/NC_010338_vs_NC_014100.fastani --fragLen $fragLen -k $kmerSize --minFraction $minFrac


fastANI -q $input_dir/NC_011916.fna -r $input_dir/NC_002696.fna -o $fastani_dir/NC_011916_vs_NC_002696.fastani --fragLen $fragLen -k $kmerSize --minFraction $minFrac


fastANI -q $input_dir/NC_011916.fna -r $input_dir/NC_010338.fna -o $fastani_dir/NC_011916_vs_NC_010338.fastani --fragLen $fragLen -k $kmerSize --minFraction $minFrac


fastANI -q $input_dir/NC_011916.fna -r $input_dir/NC_014100.fna -o $fastani_dir/NC_011916_vs_NC_014100.fastani --fragLen $fragLen -k $kmerSize --minFraction $minFrac


fastANI -q $input_dir/NC_014100.fna -r $input_dir/NC_002696.fna -o $fastani_dir/NC_014100_vs_NC_002696.fastani --fragLen $fragLen -k $kmerSize --minFraction $minFrac


fastANI -q $input_dir/NC_014100.fna -r $input_dir/NC_010338.fna -o $fastani_dir/NC_014100_vs_NC_010338.fastani --fragLen $fragLen -k $kmerSize --minFraction $minFrac


fastANI -q $input_dir/NC_014100.fna -r $input_dir/NC_011916.fna -o $fastani_dir/NC_014100_vs_NC_011916.fastani --fragLen $fragLen -k $kmerSize --minFraction $minFrac
