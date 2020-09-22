
##----------------------------------------------------------------
##                           1. Setup                           --
##----------------------------------------------------------------


## some default args that will be overwritten in the next line, but here for debugging
args <-
  c(
    "data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat", 
    20,
    "forreal",
    0
  )

#data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat
#data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat
#data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max5samp.binary.non_synonymous.dat )


## validation means use of Sanger individuals

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


## define whether a gene is being analysed
counts_logic <- FALSE
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
    )
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




## harmonise names

data_pathogen$outcome <- str_replace_all(data_pathogen$outcome, ":|\\*|\\?", "_") 

names(data_pathogen$data) <- str_replace_all(names(data_pathogen$data), ":|\\*|\\?", "_")


## define outcome
outcome <- data_pathogen$outcome

## define in plan which outcomes will have a G2G analysis
sm <-
  summarise_pathogen(data = data_pathogen,
                     format = "df",
                     filename = files_coverage)

PATHOGEN_PREDICTOR <- sm %>% 
  filter(value >= outcome_thresh & value <= (1 - outcome_thresh)) %>% 
  pull(key)
    


## make sure only phenotypes needed are in data_pathogen
data_pathogen_df <- data_pathogen$data[, c("id", PATHOGEN_PREDICTOR)]


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





##----------------------------------------------------------------
##                       5. Covar Data                          --
##----------------------------------------------------------------
source(glue::glue("{DIR_SRC}/function-covar.R"))

## WHEN DEFAULT

## this file still has all the PCs in from older data. do not use!
file_clinical <-
  glue::glue("{DIR_COVAR}/EBVG2G.covar")
data_covar_raw <- load_prepare_covar(file = file_clinical)

## separate into discrete and numeric vocariates
data_covar <-
  data_covar_raw %>% select(id, sex, age, rna)

## add info about type 1 and 2
## The values lie between -1 (type2) and 1 (type1) and have been calculated for EBNA-2 and EBNA-3s.
tmp_type <-
  fread(glue::glue("{DIR_PATHOGEN}/../clon_adj_tot_len_reads.csv")) %>% select(Sample, tot) %>% rename(id = Sample, type_1_2 = tot)

## add to covar_numeric
data_covar <- data_covar %>% left_join(tmp_type)


##----------------------------------------------------------------
##                       6. Pathogen Data                       --
##                        (continuation)                        --
##----------------------------------------------------------------



## calculate PCs or pPCs

PC <- pathogen_pca(
  path_tree = NULL,
  data_trait = data_pathogen_grm$data[, c("id", data_pathogen_grm$outcome)],
  ppca = FALSE,
  pca = TRUE,
  n.pc = n.pc,
  debugging = debugging,
  validation = validation
)


## add PCs to numeric covars
data_covar <-
  data_covar %>% left_join(PC, by = c("id" = "id"))




##----------------------------------------------------------------
##                       6. combine data                        --
##----------------------------------------------------------------


data_pathogen_covars <- data_pathogen_df %>% 
  left_join(data_covar, by = c("id" = "id")) %>%
  mutate(rna_log = log(rna))

##----------------------------------------------------------------
##                       6. run pGWAS                           --
##----------------------------------------------------------------

COVARS <- c( "type_1_2", "PC1", "PC2", "PC3", 
             "PC4", "PC5", "PC6")
COVARS_collapse <- paste(COVARS, collapse = "+")

OUTCOME <-  "rna_log"



## all coavariates
raw_model <- purrr::map_dfr(PATHOGEN_PREDICTOR, function(.x) {
  
  dat_tmp <- na.omit(data_pathogen_covars[, c(OUTCOME, .x, COVARS)])
  form <- as.formula(glue::glue("{OUTCOME} ~ {.x} + {COVARS_collapse}"))
  (lm(form , data = dat_tmp)) %>% broom::tidy() %>% mutate(pathogen = .x)
  
})

aa_model <- raw_model %>% filter(pathogen == term) %>% arrange(p.value)



## store out -----------------------------------------------------------

write_csv(aa_model, path = glue::glue("{DIR_SCRATCH}/pGWAS_viralload_{NAM}.csv"))





## only results ------------------------------

if (FALSE){
  
  

  
  PATHOGEN_PREDICTOR_SELECT <- c("BDLF1", "BCRF1", "BALF5", "BLLF3", "BORF1", "BBRF1", "BDLF4")
  PATHOGEN_PREDICTOR_SELECT <- c("BRLF1_p.Glu377Ala")
  
  raw_model <- purrr::map_dfr(PATHOGEN_PREDICTOR_SELECT, function(.x) {
    
    dat_tmp <- na.omit(data_pathogen_covars[, c(OUTCOME, .x, COVARS)])
    form <- as.formula(glue::glue("{OUTCOME} ~ {.x} + {COVARS_collapse}"))
    (lm(form , data = dat_tmp)) %>% broom::tidy() %>% mutate(pathogen = .x)
    
  })
  
  aa_model <- raw_model %>% filter(pathogen == term) %>% arrange(p.value)
  
}



