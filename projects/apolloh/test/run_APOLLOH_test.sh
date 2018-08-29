#!/bin/bash

### CHANGE TO ROOT INSTALL DIR ###
cd ../

### SETUP VARIABLES ###
installDir=$PWD/
outDir=test/
inputFile=test/infile.txt
copyNumberFile=test/cnfile.txt

### SETUP PATH ###
export LD_LIBRARY_PATH=$installDir/MATLAB_Component_Runtime/v77/runtime/glnxa64/:$installDir/MATLAB_Component_Runtime/v77/sys/os/glnxa64/:$installDir/MATLAB_Component_Runtime/v77/bin/glnxa64/:$LD_LIBRARY_PATH


##### EXECUTE APOLLOH ####
#Arguments:
#1) Tab-delimited input file containing allelic counts from the tumour at positions determined as heterozygous from the normal genome
#2) Tab-delimited input copy number segment prior file. The accepted format is the output from HMMcopy, a read-depth for analyzing copy number in tumour-normal sequenced genomes.
#3) Parameter intialization file is a matlab binary (.mat) file. This file contains model and setting paramters necessary to run the program. There are specific parameters for Illumina and SOLiD sequencing data. This example is for SOLiD data.
#4) Tab-delimited output file for position-level results.
#5) Tab-delimited output file storing converged parameters after model training using Expectation Maximization (EM) algorithm.
time $installDir/bin/apolloh $inputFile $copyNumberFile $installDir/parameters/apolloh_K18_params_SOLiDskew_stromalRatio_Hyper10k.mat $outDir/test_params.txt $outDir/test_loh.txt

### Sort the output file by chromosome and position ###
sort -k1,1n -k2,2n $outDir/test_loh.txt > $outDir/test_loh.txt.tmp; mv $outDir/test_loh.txt.tmp $outDir/test_loh.txt;

### Create a segment file from the output file using a perl script ###
$installDir/scripts/createSingleSegFileFromAPOLLOH.pl -i $outDir/test_loh.txt -o $outDir/test_segs.txt -calls 1;
