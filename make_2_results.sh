#!/bin/bash

## ///////////////
## PARAMETERS
## ///////////////

MODE="forreal"      ## "debugging", "forreal", "random"
CORES=25              ## number of cores
TAXID="t1"            ## there are two taxids: t1 and t2.

## how many PATHOGEN is determined in setup.R (search for object NAMS)

PATHOGEN=(data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat
          data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat )

OUTCOMEMODEL=( binary binary)


## run to get more stricter QC SNPs
## run only for 1 pathogen
rm -rf data/1_processed/*
Rscript --vanilla src/post_QC.R  data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat $CORES $MODE 0.1
rm -rf data/1_processed/*



## ///////////////////////////
## run over PATHOGEN
## ///////////////////////////

for index in ${!PATHOGEN[*]} 
do

echo "Started phenotype $index"

printf "   %s\n" $index

## run results -------------------------------------------
Rscript --vanilla src/results.R  ${PATHOGEN[$index]} $TAXID $MODE "gcta"

echo "Done with $index phenotype"

## clean up data/1_processed/
rm -rf data/1_processed/*

done
