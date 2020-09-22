## produce a more lenient QC for the host genotype data
## --------------------------------------------------

##----------------------------------------------------------------
##                           1. Setup                           --
##----------------------------------------------------------------


## some default args that will be overwritten in the next line, but here for debugging
args <-
  c(
    "data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat", 
    20,
    "forreal",
    0.1
  )

args <- commandArgs(trailingOnly = TRUE)

PATHOGEN <- args[1]
n_cores <- as.numeric(args[2]) ## number of threads used
debugging <- FALSE
random_outcome <- FALSE
outcome_thresh <- 0.1

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


## make callrate thresholds more strict

callrate.snps.thresh <- 0.01
callrate.inds.thresh <- 0.03

## QC:
HOST_QC <-
  glue::glue("{DIR_PROCESSED}/{NAM_HOST}_QC")

## select individuals first
## run parts of src/run-gcta.R first
system(
  glue::glue(
    "{PLINK} --bfile {HOST} --make-bed --out {HOST_QC}_intermediate1 --keep {DIR_PROCESSED}/{str_forgcta}_idlist.txt"
  )
)


## For GCTA MLMA --------------
## apply callrate snps
system(
  glue::glue(
    "{PLINK} --bfile {HOST} --geno {callrate.snps.thresh} --maf {maf.thresh} --hwe {hwe.thresh} --make-bed --out {HOST_QC}_intermediate2"
  )
)
## apply callrate individuals
system(
  glue::glue(
    "{PLINK} --bfile {HOST_QC}_intermediate2 --mind {callrate.inds.thresh} --maf {maf.thresh} --hwe {hwe.thresh} --make-bed --out {HOST_QC}"
  )
)



## remove fam + bed
## QC:
## copy bim to DIR_SCRATCH
## remove all DIR_PROCESSED

fs::file_copy(glue::glue("{HOST_QC}.bim"), glue::glue("{DIR_SCRATCH}/{NAM_HOST}_QC.bim"), overwrite = TRUE)
