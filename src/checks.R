


############################################################################
############################################################################
###                                                                      ###
###                             SANITYCHECKS                             ###
###                                                                      ###
############################################################################
############################################################################

## Perform some sanity checks on data/2_results
## 0. Setup
## 1. Does file exists that should exist?
## 2. Does file exists that should not exist?    
## 3. if file listed, does it actually exist?




##////////////////////////////////////////////////////////////////
##                           0. Setup                           //
##////////////////////////////////////////////////////////////////

#PATHOGEN=(/home/rueger/G2G-EBV/data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max5samp.binary.non_synonymous.dat 
#          /home/rueger/G2G-EBV/data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat
#          /home/rueger/G2G-EBV/data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_aa_variant_matrix.binary.dat)

## some default args that will be overwritten in the next line, but here for debugging
args <-
  c(
    "data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat",
    "t1",
    "forreal"
  )


args <- commandArgs(trailingOnly = TRUE)


## some default settings
PATHOGEN <- args[1]
taxid <- args[2]
debugging <-
  ifelse(args[3] == "debugging", TRUE, FALSE) ## if "debugging" make datasets small for checking if it runs
random_outcome <-
  ifelse(args[3] == "random", TRUE, FALSE) ## if "random", create random outcome vectors with differing frequency of 1's and 0's

## This R script creates all the variables
## takes PATHOGEN as input!
DIR_SRC <- here::here("src")
source(glue::glue("{DIR_SRC}/setup.R"))


## add string to NAM, to indicate debugging, random mode
if (debugging)
  NAM <- glue::glue("{NAM}_debugging")
if (random_outcome)
  NAM <- glue::glue("{NAM}_random")


## read in plan
plan <-
  read_delim(glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"), delim = " ") %>% mutate(str.out.gcta = str_replace(str.out.gcta, "/work/backup/gr-fe/rueger/G2G-EBV/data/2_results/", glue::glue("{DIR_SCRATCH}/")))


##////////////////////////////////////////////////////////////////
##            CHECK 1: file exists that should exist            //
##////////////////////////////////////////////////////////////////
out <- data.frame()

check.NA.exists <-
  plan %>% filter(is.na(str.out.gcta) &
                    include.gcta) %>% select(global.outcome, outcome, str.out.gcta)

if (nrow(check.NA.exists) > 0){
  message1 <- "Some outcome did not run"
  ## stop(message)
  
  check.NA.exists <- data.frame(check.NA.exists,error =  message1)
  out <- rbind(out, check.NA.exists)
}


##////////////////////////////////////////////////////////////////
##            CHECK 2: file exists that should not exist        //
##////////////////////////////////////////////////////////////////

files <-
  fs::dir_ls(DIR_SCRATCH) %>%
  str_subset(NAM) %>%
  str_ignore("subset.txt") %>%
  str_ignore("validation") %>%
  str_subset("gcta") %>%
  str_ignore("clumped") %>%
  str_ignore("log") %>%
  str_ignore("nosex")
  


files_should_not_compute <- files[!(files %in% plan$str.out.gcta)]

out <- data.frame()


if (length(files_should_not_compute) > 0){
  message2 <- "Outcome computed that should not compute"
  ## stop(message)
  
  check.filesnot.exists <- data.frame(files_should_not_compute, error =  message2)
  out <- rbind(out, check.filesnot.exists)
}


##////////////////////////////////////////////////////////////////
##     CHECK 3: if file listed, does it actually exist?         //
##////////////////////////////////////////////////////////////////

check.file.exists <-
  plan %>% filter(!is.na(str.out.gcta)) %>% filter(!file.exists(str.out.gcta)) %>% select(global.outcome, outcome, str.out.gcta)

if (nrow(check.file.exists) > 0){
  message3 <- "This outcome did not compute (no file exists)"
  ## stop(message2)
  check.file.exists <- data.frame(check.file.exists, error = message3)
  out <- rbind(out, check.file.exists)
  
}

write_delim(out, path = glue::glue("{DIR_OUTPUT}/checks_{NAM}.log"), delim = " " )
