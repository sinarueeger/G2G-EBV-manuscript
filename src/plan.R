


###########################################################################
###########################################################################
###                                                                     ###
###                              MAKE PLAN                              ###
###                                                                     ###
###########################################################################
###########################################################################
## this script lays out the plan of what is analysed
###########################################################################


#make  df with
#[outcome, outcome.freq, predictors, model, str.out.gcta, str.out.lmm]

##----------------------------------------------------------------
##                           1. Setup                           --
##----------------------------------------------------------------

## some default args that will be overwritten in the next line, but here for debugging
args <-
  c(
    "data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat", 
    "binary",
    "lmm",
    "t1",
    "forreal",
    ""
  )

args <- commandArgs(trailingOnly = TRUE)

PATHOGEN <- args[1]
model_outcome <- args[2]
model_method <- args[3]
taxid <- args[4]
debugging <-
  ifelse(args[5] == "debugging", TRUE, FALSE) ## if "debugging" make datasets small for checking if it runs
random_outcome <-
  ifelse(args[5] == "random", TRUE, FALSE) ## if "random", create random outcome vectors with differing frequency of 1's and 0's
n_random <- 20 ## number of random vecs (i.e. number of GWASs)
gene_str <- args[6]
gene_test <-
  ifelse(gene_str != "all", TRUE, FALSE) ## if string is passed on and not ""

## This R script creates all the variables
DIR_SRC <- here::here("src")
source(glue::glue("{DIR_SRC}/setup.R"))
## Here is a list of the variables it creates that we are recycling later
## - HOST: host file string, with path
## - NAM: thats the PATHOGEN file string, but without path
## - PATHOGEN: same as NAM but no path
## - DIR_SRC: where this file is stored
## - DIR_SCRATCH: where all results are stored

## rename it before changing it
NAM_raw <- NAM

## add string to debugging
if (debugging) {
  NAM <- glue::glue("{NAM}_debugging")
}

## add string to random
if (random_outcome) {
  NAM <- glue::glue("{NAM}_random")
}


##----------------------------------------------------------------
##                      2. PATHOGEN DATA                        --
##----------------------------------------------------------------


if (random_outcome)
{
  ## when random pathogen vector
  ## generate random vectors with frequency between 0 and 1
  sm <-
    data.frame(
      key = paste0("random_vec_", 1:n_random),
      value = seq(0, 1, length.out = n_random + 1)[-1]
    )
  
} else {
  ## get functions to treat pathogen
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
  
  ## files_coverage provides us with a quality measure of the pathogen files.
  
  ## load raw files ---
  data_pathogen_raw <-
    load_prepare_pathogen(
      files_coverage = files_coverage,
      files_data = PATHOGEN,
      debugging = debugging,
      path_gene_length = FILE_PATH_COUNTS,
      counts2frac = model_outcome == "continuous" & sum(str_detect(PATHOGEN, c("gene", "counts"))) == 2 ## in case it is a gene and counts logic, gene or not
    ) 
  
  outcome <- data_pathogen_raw$outcome
  
  sm <-
    summarise_pathogen(data = data_pathogen_raw,
                       format = "df",
                       filename = files_coverage)
  
}


##----------------------------------------------------------------
##                       3. Create Plan                         --
##----------------------------------------------------------------

if (gene_test)
{
  ## if specific gene string should be tested, create a subset and then filter sm
  outcome_subset <- outcome %>% str_subset(gene_str)
  sm <- sm %>% filter(key %in% outcome_subset)
}

## Assemble plan
plan <-
  data.frame(
    global.outcome = NAM,
    outcome = sm$key,
    outcome.freq = sm$value,
    taxid = taxid,
    model = model_method,
    model_outcome = model_outcome,
    covars.gcta = NA,
    str.out.gcta = NA,
    covars.null = NA,
    str.out.null = NA,
    host = HOST
  )

if (debugging)
{
  ## if debugging, just select first three
  plan <- plan[1:3, ]
}

## write out plan
write_delim(plan,
            path = glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"),
            delim = " ")
