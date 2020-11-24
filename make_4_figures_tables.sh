#!/bin/bash

## ///////////////
## PARAMETERS
## ///////////////

MODE=("forreal"        ## "debugging", "forreal", "random"
       "forreal"
       "forreal"
       "forreal"
      )
CORES=1                ## number of cores (nproc --all) ## cause running GWAS at the same time
MODEL="lmm"             ## meaning - mixed model. there is no other option implemented.
TAXID="t1"              ## there are two taxids: t1 and t2.

## how many PATHOGEN is determined in setup.R (search for object NAMS)
PATHOGEN=(/home/rueger/G2G-EBV/data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat
          /home/rueger/G2G-EBV/data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat
          )

## THRESHOLD outcome
THRESH=( 0.1 0 ) ## outcome frequency threshold, only select outcomes with freq >= THRESH

OUTCOMEMODEL=( binary binary ) ## binary

DIR_MS=data/5_manuscript

mkdir -p $DIR_MS


## ///////////////////////////
## run descriptive
## ///////////////////////////

for index in ${!PATHOGEN[*]} 
do

      echo "Started phenotype $index"
      printf "%s\n" $index

      ## mimics G2G data prep & creates all data in 4_descript
      Rscript --vanilla report/descript-generate-data.R  ${PATHOGEN[$index]} $MODEL ${OUTCOMEMODEL[$index]} $TAXID $CORES ${MODE[$index]} ${THRESH[$index]}

      ## input: 4_descript, output 4_descript (but higher level)
      Rscript --vanilla report/descript-cohort.R  ${PATHOGEN[$index]} ## stores out to 4_descript

      echo "Done with $index phenotype"

done


## ///////////////////////////
## run all figures + tables
## ///////////////////////////
## input: 4_descript, output 5_manuscript

Rscript --vanilla report/figures-tables.R
