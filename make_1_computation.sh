#!/bin/bash

## ///////////////
## PARAMETERS
## ///////////////

MODE="forreal"          ## "debugging", "forreal", "random"
CORES=25                ## when on SLURM $SLURM_NTASKS     ## number of cores (nproc --all)
MODEL="lmm"             ## meaning - mixed model. there is no other option implemented.
TAXID="t1"              ## there are two taxids: t1 and t2.

## how many PATHOGEN is determined in setup.R (search for object NAMS)
PATHOGEN=(data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat 
          data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat )

OUTCOMEMODEL=( binary binary ) ## binary

## THRESHOLD outcome
THRESH=( 0.1 0.04 ) ## outcome frequency threshold, only select outcomes with freq >= THRESH

## whether a specific outcome should be computed
TESTGENE=all ## "all" or some regex expression, e.g. BBLF2 or ## outcome <- c("EBNA_3B_t1_alt:p.Arg338Arg", "BBLF2_BBLF3:p.Val204Ile")


## ///////////////////////////
## run over PATHOGEN
## ///////////////////////////

for index in ${!PATHOGEN[*]} 
do

echo "Started phenotype $index"
printf "%s\n" $index

## clean up data/1_processed/
rm -rf data/1_processed/*

## execute plan -------------------------------------
#if [ "$MODE" == "debugging" ]; then

echo "STARTED PLAN //////////// "
Rscript --vanilla src/plan.R ${PATHOGEN[$index]} ${OUTCOMEMODEL[$index]} $MODEL $TAXID $MODE $HOST $TESTGENE
echo "//////////// FINISHED PLAN "
echo $TESTGENE
#fi

## run gcta -----------------------------------------
Rscript --vanilla src/make-gcta.R  ${PATHOGEN[$index]} $CORES $MODE ${THRESH[$index]}
rm -rf data/1_processed/*

## run pGWAS ---------------------------------------
Rscript --vanilla src/make-viralload-pGWAS.R  ${PATHOGEN[$index]} $CORES $MODE ${THRESH[$index]}
rm -rf data/1_processed/*

echo "Done with $index phenotype"

done
