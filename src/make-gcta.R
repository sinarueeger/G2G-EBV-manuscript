############################################################################
############################################################################
###                                                                      ###
###                           G2G-EBV ANALYSIS                           ###
###                              USING GCTA                              ###
###                                                                      ###
############################################################################
############################################################################


## This script aims to run an (almost) identical analysis as in make-gaston.R, but
## with GCTA instead of the R-package gaston
## ############################################################################


##::::::::::::::::::::::::::::::::
## The input are several variables defined in
## src/setup.R + the variable args passed on from the bash file.
## The output are png's and txt's files in the DIR_SCRATCH folder
##::::::::::::::::::::::::::::::::


## 1. Setup
## 2. Pathogen data
## 3. Subset data for individuals         
## 4. Host data
## 5. Covar data
## 6. Pathogen data (continuation)
## 7. GRM
## 8. check order of files        
## 9. G2G


##----------------------------------------------------------------
##                           1. Setup                           --
##----------------------------------------------------------------


## some default args that will be overwritten in the next line, but here for debugging
args <-
  c(
    "data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat", 
    20,
    "forreal",
    0.04
  )

args <- commandArgs(trailingOnly = TRUE)

PATHOGEN <- args[1]
n_cores <- as.numeric(args[2]) ## number of threads used
debugging <-
  ifelse(args[3] == "debugging", TRUE, FALSE) ## if "debugging" make datasets small for checking if it run
random_outcome <-
  ifelse(args[3] == "random", TRUE, FALSE) ## if "random", create random outcome vectors with differing frequency of 1's and 0's
outcome_thresh <-
  as.numeric(args[4]) ## threshold for outcome frequency (should be 0.3 for AA and 0.1 for genes)


## Set paths ---

## This R script creates all the variables
DIR_SRC <- here::here("src")
source(glue::glue("{DIR_SRC}/setup.R"))

## when running cohort-description.R
## NAM2 <- str_replace(PATHOGEN, "/home/rueger/G2G-EBV/data/0_raw/pathogen/", "") %>% str_replace("/", "_")
## DIR_PROCESSED <- glue::glue("/home/rueger/G2G-EBV/data/4_descript/{NAM2}")
## fs::dir_create(DIR_PROCESSED)

## add string to debugging
if (debugging)
{
  NAM <- glue::glue("{NAM}_debugging")
  
  ## cause only 3 outcomes preselected, drop outcome_thresh
  outcome_thresh <- 0
}

## Read in the plan generated in DIR_SRC/plan.R
if (file.exists(glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt")))
{
  plan <-
    read_delim(glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"), delim = " ")
} else{
  stop("no plan for this phenotype")
}


model_outcome <- unique(plan$model_outcome)
model_method <- unique(plan$model)
taxid <- unique(plan$taxid)

stopifnot(length(c(unique(plan$model_outcome), unique(plan$model), unique(plan$taxid))) == 3)


## DEFINE HOST
f_check_length <- function(x)
  if (length(x) > 1)
    return(FALSE)

## Define HOST data
HOST <-
  plan %>% select(host) %>% distinct() %>% assertr::verify(nrow(.) == 1) %>% pull(host)
NAM_HOST <- fs::path_file(HOST)##  glue::glue("{(DIR_HOST)}/"))

## needed to name all files
str_forgcta <- "forgcta"
if (debugging)
  str_forgcta <- "forgcta_debugging"
#if(random_outcome) NAM <- glue::glue("{NAM}_random")

## define whether a gene is being analysed
counts_logic <- model_outcome == "continuous" & sum(str_detect(PATHOGEN, c("gene", "counts"))) == 2
## in case it is a gene and counts logic, gene or not

n.pc <- 6

##----------------------------------------------------------------
##                       2. Pathogen Data                       --
##----------------------------------------------------------------

source(glue::glue("{DIR_SRC}/function-pathogen.R"))

## define coverage filenames ---
files_coverage <- case_when(
  str_detect(PATHOGEN, "t1") ~ as.character(
    glue::glue("{DIR_PATHOGEN}/../NC_009334_with_t1_alts.covstats.dat")
  ),
  str_detect(PATHOGEN, "t2") ~ as.character(
    glue::glue("{DIR_PATHOGEN}/../NC_007605_with_t2_alts.covstats.dat")
  ),
  TRUE ~ NA_character_
)

## load raw files ---
data_pathogen_raw <-
  load_prepare_pathogen(
    files_coverage = files_coverage,
    files_data = PATHOGEN,
    debugging = debugging,
    path_gene_length = FILE_PATH_COUNTS,
    ## in case its gene
    counts2frac = counts_logic ## in case it is a gene and counts logic, gene or not
  ) ## logic, gene or not


## apply QC to pathogen data ---
data_pathogen <- apply_qc_pathogen(data_pathogen_raw)
## data_pathogen is a list!
## data_pathogen$data
## data_pathogen$outcome


## define outcome
outcome <- data_pathogen$outcome


## define in plan which outcomes will have a G2G analysis
if (counts_logic)
{
  ## filter for outcome_thresh first
  plan <- plan %>% mutate(include.gcta = TRUE)
  
} else{
  ## filter for outcome_thresh first
  plan <-
    plan %>% mutate(include.gcta = (
      outcome.freq >= outcome_thresh &
        outcome.freq <= (1 - outcome_thresh)
    ))
  
}



## we have to take the intersect between the QC done before for outcome, and the QC based on the thresholding for frequencies)
outcome <- intersect(outcome, plan$outcome[plan$include.gcta])


## make sure only phenotypes needed are in data_pathogen
data_pathogen$data <- data_pathogen$data[, c("id", outcome)]


## needed for later, to loop through the outcomes
n.outcome <- length(outcome)


## calculate PCs from pathogen variants
data_pathogen_grm <-
  apply_qc_pathogen(
    load_prepare_pathogen(
      files_coverage = files_coverage,
      files_data = PATHOGEN_FORGRM,
      debugging = debugging,
      counts2frac = FALSE
    )
  )





##---------------------------------------------------------------
##                 3. Subset data for individuals              --
##---------------------------------------------------------------
## ... for individuals that have pathogen data available

## create a id list, so that genotype data can be extracted
idlist <-
  data.frame(id = data_pathogen$data$id, id2 = data_pathogen$data$id)
write_delim(
  idlist,
  path = glue::glue("{DIR_PROCESSED}/{str_forgcta}_idlist.txt"),
  col_names = FALSE,
  delim = " "
)


##----------------------------------------------------------------
##                         4. Host Data                         --
##----------------------------------------------------------------
## see also https://onlinelibrary.wiley.com/doi/epdf/10.1002/mpr.1608
## {HOST}

## QC:
HOST_QC <-
  glue::glue("{DIR_PROCESSED}/{NAM_HOST}_QC")

if (debugging) {
  ## when DEBUGGING
  
  HOST_QC <- glue::glue("{HOST_QC}_debugging")
  
  ## select individuals first
  system(glue::glue(
    "{PLINK} --bfile {HOST} --make-bed --out {HOST_QC}_intermediate1 --keep {DIR_PROCESSED}/{str_forgcta}_idlist.txt"
  ))
  
  ## apply callrate snps
  system(glue::glue(
    "{PLINK} --bfile {HOST_QC}_intermediate1 --geno {callrate.snps.thresh} --from-bp {from_to_debugging[1]} --to-bp {from_to_debugging[2]} --chr {chr_debugging} --maf {maf.thresh} --hwe {hwe.thresh} --make-bed --out {HOST_QC}_intermediate2"
  ))
  
  ## apply callrate individuals
  system(glue::glue(
    "{PLINK} --bfile {HOST_QC}_intermediate2 --from-bp {from_to_debugging[1]} --to-bp {from_to_debugging[2]} --chr {chr_debugging} --mind {callrate.inds.thresh} --maf {maf.thresh} --hwe {hwe.thresh} --make-bed --out {HOST_QC}"
  ))
  
  
} else {
  ## by DEFAULT
  
  ## select individuals first
  system(
    glue::glue(
      "{PLINK} --bfile {HOST} --make-bed --out {HOST_QC}_intermediate1 --keep {DIR_PROCESSED}/{str_forgcta}_idlist.txt"
    )
  )
  
  
  ## For GCTA MLMA --------------
  ## apply callrate snps
  system(glue::glue(
    "{PLINK} --bfile {HOST_QC}_intermediate1 --geno {callrate.snps.thresh} --maf {maf.thresh} --hwe {hwe.thresh} --make-bed --out {HOST_QC}_intermediate2"
  ))
  ## apply callrate individuals
  system(glue::glue(
    "{PLINK} --bfile {HOST_QC}_intermediate2 --mind {callrate.inds.thresh} --maf {maf.thresh} --hwe {hwe.thresh} --make-bed --out {HOST_QC}"
  ))
  
  ## For GRM -------------------
  system(glue::glue(
    "{PLINK} --bfile {HOST_QC}_intermediate1 --geno 0.01  --maf {maf.thresh} --hwe {hwe.thresh} --make-bed --out {HOST_QC}_intermediate2_GRM"
  ))
  ## apply callrate individuals
  system(glue::glue(
    "{PLINK} --bfile {HOST_QC}_intermediate2_GRM --mind 0.03 --maf {maf.thresh} --hwe {hwe.thresh} --make-bed --out {HOST_QC}_GRM"
  ))
  
}


## read the fam file (for later purposes)
fam.file <-
  read_delim(glue::glue("{HOST_QC}.fam"),
             delim = "\t",
             col_names = FALSE)


##----------------------------------------------------------------
##                       5. Covar Data                          --
##----------------------------------------------------------------
source(glue::glue("{DIR_SRC}/function-covar.R"))

## WHEN DEFAULT

## this file still has all the PCs in from older data. do not use!
file_clinical <-
  glue::glue("{DIR_COVAR}/EBVG2G.covar")
data_covar <- load_prepare_covar(file = file_clinical)

## add info from 17819 and 31570!!!
#  Sina asked: IDs `17819` and `31570` of the SHCS individuals that have missing covariates.
#  Jacques answered: Strangely enough, both are born in 1974, and both are male.
data_covar_additional <-
  tibble(
    id = c(17819, 31570),
    sex = c(1, 1),
    yearofbirth = as.Date(c("1974-01-01", "1974-01-01"))
  ) %>% mutate(age = as.period(lubridate::interval(yearofbirth, "2004-09-23"), unit = "year")$year)
## https://gist.github.com/mmparker/7254445

#> data_covar[sample(1:nrow(data_covar), 3), c(1, 6, 7)]
## A tibble: 3 x 3
# fid   sex   age
# <int> <int> <dbl>
# 1 31627     1  32.8
# 2 25839     2  37.2
# 3 16273     1  26.7
##61627: 1972 --- 25839: 1962 --- 16273: 1969
## Age corresponds to the time of consent and registration into the cohort (23 sep 2004)


## separate into discrete and numeric vocariates
covar_discrete <-
  data_covar %>% select(id, sex) %>% bind_rows(data_covar_additional  %>% select(id, sex))
covar_numeric <-
  data_covar %>% select(id, age) %>% bind_rows(data_covar_additional  %>% select(id, age))

## add info about type 1 and 2
## The values lie between -1 (type2) and 1 (type1) and have been calculated for EBNA-2 and EBNA-3s.
tmp_type <-
  fread(glue::glue("{DIR_PATHOGEN}/../clon_adj_tot_len_reads.csv")) %>% select(Sample, tot) %>% rename(id = Sample, type_1_2 = tot)

## add to covar_numeric
covar_numeric <- covar_numeric %>% left_join(tmp_type)


## turn into same order than genotyped data
covar_numeric <-
  left_join(fam.file %>% select(X1, X2),
            covar_numeric,
            by = c("X1" = "id"))

covar_discrete <-
  left_join(fam.file %>% select(X1, X2),
            covar_discrete,
            by = c("X1" = "id"))

## sanity check
stopifnot(nrow(covar_discrete) == 268)
stopifnot(nrow(covar_numeric) == 268)





##----------------------------------------------------------------
##                       6. Pathogen Data                       --
##                        (continuation)                        --
##----------------------------------------------------------------

## Create pathogen DATA FRAME (before it was a list)
data_pathogen_df <-
  left_join(fam.file %>% select(X1, X2), data_pathogen$data, by = c("X1" = "id"))


## phen file
## Format ID, FAMID, Pheno
write_delim(
  data_pathogen_df,
  path = glue::glue("{DIR_PROCESSED}/{str_forgcta}_{NAM}.pheno"),
  col_names = FALSE,
  delim = " "
)


## Sanitycheck
## check if identical columns
ped <- read_tsv(glue::glue("{HOST_QC}.fam"), col_names = FALSE)
stopifnot(identical(as.character(ped$X1), as.character(data_pathogen_df$X1)))


## calculate PCs or pPCs
##path_tree_EBNA_3A <- glue::glue("{DIR_PATHOGEN}/../EBNA_3A.hiv.msa.manually_trimmed.fasta.treefile") # ape::read.tree(path_tree_EBNA_3A)
##path_tree_EBNA_2 <- glue::glue("{DIR_PATHOGEN}/../EBNA_2.hiv.msa.manually_trimmed.fasta.treefile")

PC <- pathogen_pca(
  path_tree = NULL,
  data_trait = data_pathogen_grm$data[, c("id", data_pathogen_grm$outcome)],
  ppca = FALSE,
  pca = TRUE,
  n.pc = n.pc,
  debugging = debugging
)


## add PCs to numeric covars
covar_numeric <-
  covar_numeric %>% left_join(PC, by = c("X1" = "id"))


## write out
write_delim(
  covar_discrete,
  glue::glue("{DIR_PROCESSED}/{str_forgcta}_covar_discrete.covar"),
  delim = " ",
  col_names = FALSE
)
write_delim(
  covar_numeric,
  glue::glue("{DIR_PROCESSED}/{str_forgcta}_covar_numeric.qcovar"),
  delim = " ",
  col_names = FALSE
)


## viralload
viralload <-
  data_covar %>% select(id, rna) %>% mutate(rna_log = log(rna)) %>% select(-rna)
viralload <- left_join(fam.file %>% select(X1, X2),
                       viralload,
                       by = c("X1" = "id"))


## write out
write_delim(
  viralload,
  glue::glue("{DIR_PROCESSED}/{str_forgcta}_viralload.pheno"),
  delim = " ",
  col_names = FALSE
)

##---------------------------------------------------------------
##                           7. GRMs                           --
##---------------------------------------------------------------

system(
  glue::glue(
    "{GCTA} --bfile {HOST_QC}_GRM --autosome --maf 0.05 --make-grm --out {DIR_PROCESSED}/{str_forgcta}_GRM --thread-num {n_cores}"
  )
)

##---------------------------------------------------------------
##                  8. check order of files                    --
##---------------------------------------------------------------


host_id <- read_tsv(glue::glue("{HOST_QC}.fam"), col_names = FALSE)
pheno_id <-
  read_delim(
    glue::glue("{DIR_PROCESSED}/{str_forgcta}_{NAM}.pheno"),
    col_names = FALSE,
    delim = " "
  )
covar_id <-
  read_delim(
    glue::glue("{DIR_PROCESSED}/{str_forgcta}_covar_discrete.covar"),
    col_names = FALSE,
    delim = " "
  )
qcovar_id <-
  read_delim(
    glue::glue("{DIR_PROCESSED}/{str_forgcta}_covar_numeric.qcovar"),
    col_names = FALSE,
    delim = " "
  )


## Sanitychecks
stopifnot(identical(as.character(host_id$X1), as.character(pheno_id$X1)))
stopifnot(identical(as.character(host_id$X1), as.character(covar_id$X1)))
stopifnot(identical(as.character(host_id$X1), as.character(qcovar_id$X1)))
stopifnot(identical(as.character(host_id$X1), as.character(viralload$X1)))


##----------------------------------------------------------------
##                           9. GCTA                            --
##----------------------------------------------------------------


## Sanitychecks
## the order of outcome here is important!
stopifnot(identical(names(data_pathogen_df)[-c(1:2)], outcome))


## Loop through all outcomes

for (counter_within_outcome in 1:n.outcome) {
  ## string to write out
  str_out <-
    glue::glue("{model_method}_gcta_{NAM}_{outcome[counter_within_outcome]}")
  
  
  ## update plan
  plan[plan$outcome == outcome[counter_within_outcome], "str.out.gcta"] <-
    glue::glue("{DIR_SCRATCH}/{str_out}.mlma")
  plan[plan$outcome == outcome[counter_within_outcome], "covars.gcta"] <-
    paste(c(names(covar_discrete)[-c(1:2)], names(covar_numeric)[-c(1:2)]), collapse = " + ")
  write_delim(plan,
              path = glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"),
              delim = " ")
  
  
  
  ## Run G2G analysis
  ## removed --mlma-no-preadj-covar (If this option is specified, the covariates will be fitted together with the SNP for association test. However, this will significantly reduce computational efficiency.)
  ## info http://gcta.freeforums.net/thread/247/greml-estimating-variance-explained-snps
  
  system(
    glue::glue(
      "{GCTA} --mlma --bfile {HOST_QC} --grm {DIR_PROCESSED}/{str_forgcta}_GRM --pheno {DIR_PROCESSED}/{str_forgcta}_{NAM}.pheno --mpheno {counter_within_outcome} --out {DIR_SCRATCH}/{str_out} --thread-num {n_cores} --covar {DIR_PROCESSED}/{str_forgcta}_covar_discrete.covar --qcovar {DIR_PROCESSED}/{str_forgcta}_covar_numeric.qcovar"
    )
  )
  
}

## write out plan (this step is not really necessary)
write_delim(plan,
            path = glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"),
            delim = " ")


##----------------------------------------------------------------
##                        10. Viralload                         --
##----------------------------------------------------------------

str_out <-
  glue::glue("GWAS_viralload")

system(
  glue::glue(
    "{GCTA} --mlma --bfile {HOST_QC} --grm {DIR_PROCESSED}/{str_forgcta}_GRM --pheno {DIR_PROCESSED}/{str_forgcta}_viralload.pheno --out {DIR_SCRATCH}/{str_out} --thread-num {n_cores} --covar {DIR_PROCESSED}/{str_forgcta}_covar_discrete.covar --qcovar {DIR_PROCESSED}/{str_forgcta}_covar_numeric.qcovar"
  )
)

## out <- data.table::fread(glue::glue("{DIR_SCRATCH}/GWAS_viralload.mlma"))
