#!/bin/bash

## this runs the genomic lambda calculation. 

## ///////////////
## PARAMETERS
## ///////////////


PATHOGEN=(data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat 
          data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat
          data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat
          data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max5samp.binary.non_synonymous.dat)
          

## ///////////////
## CALC GC LAMBDA
## ///////////////

for index in ${!PATHOGEN[*]} 
do

echo "Started phenotype $index"

printf "   %s\n" $index

Rscript --vanilla src/calc-genomic-lambda.R  ${PATHOGEN[$index]}

done
