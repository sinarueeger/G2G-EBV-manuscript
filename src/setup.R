############################################################################
############################################################################
###                                                                      ###
###                                SETUP                                 ###
###                                                                      ###
############################################################################
############################################################################

## This file is needed to set up any analysis for the G2G analysis

## what it does:
## - it loads all libraries needed
## - defines all paths
## - defines extra variables

## things to pass along:
## - {taxid} ("t1" or "t2")

## #########################################################################

##////////////////////////////////////////////////////////////////
##                   setting defaults                           //
##////////////////////////////////////////////////////////////////

## setting defaults if they don't exist
if (!exists("taxid")) {
  taxid <- "t1"
  message("setting default for taxid")
}

##////////////////////////////////////////////////////////////////
##                      load libraries                          //
##////////////////////////////////////////////////////////////////

## for all files ------------------------------------------------

library(glue)
library(data.table)
library(dplyr)
library(devtools)
library(stringr)
library(readr)
library(janitor)
library(here)
library(magrittr)
library(tidyr)
library(fs)
library(purrr)
library(knitr)
library(tidyverse) ## needs to be here at the end, even though redundant to make the kable tables work..
library(lubridate)
options(datatable.fread.datatable = FALSE) ## will return from fread a dataframe and not a datatable
library(logisticPCA) ## for function logisticPCA in src/function-pathogen.R
library(missMDA) ## for function imputePCA in src/function-pathogen.R (installed 3.0-3 to get it working)
library(qqman)
library(broom)
## for make_2_results.sh ------------------------------------------
library(furrr) ## for parallelized purrr calls
library(ggman)

library(colorblindr)
# install.packages("colorspace", repos = "http://R-Forge.R-project.org")
#library(gplots) ## for heatmap, addressing with ::
#library(patchwork) using, but with :: devtools::install_github("thomasp85/patchwork")
library(ggplot2)

## get theme
devtools::source_url("https://raw.github.com/sinarueeger/ggtheme/master/theme_presentation.R")


epfl_red <- "#FF0000"
epfl_turquise <- "#00A79F"

epfl_dark_red <- "#B51F1F"
epfl_dark_turquise <- "#007480"

epfl_gray <- "#CAC7C7"
epfl_dark_gray <- "#413D3A"


library(ggGWAS) ## via github.com/sinarueeger/ggGWAS


## to sort out conflicted packages ---------------------------------

library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("chisq.test", "stats")
conflict_prefer("fisher.test", "stats")

##////////////////////////////////////////////////////////////////
##                      define paths                            //
##////////////////////////////////////////////////////////////////

## all paths are defined with here::here
## data is softlinked

## some extra paths
DIR_HOST <- here::here("data", "0_raw", "host", "geno")
DIR_PATHOGEN <- here::here("data", "0_raw", "pathogen", "results_max3_bwa_mem_hiv_18")
DIR_COVAR <- here::here("data", "0_raw", "host", "covar")
DIR_PROCESSED <-
  here::here("data", "1_processed") ## this is for all processed raw data
DIR_OUTPUT <-
  here::here("data", "3_output") ## this is where descriptive stats goes
DIR_DESCRIPT <-
  here::here("data", "4_descript") ## this is where descriptive stats goes
DIR_MS <-
  here::here("data", "5_manuscript") ## this is where ms figures and tables go

## define where results should be stored
COMPUTER <- system("uname -n", intern = TRUE)

## this is where all the written results go
DIR_SCRATCH <- case_when(
  COMPUTER %in% c("deneb1", "deneb2", "fidis") ~ here::here("data", "2_results"),
  str_detect(COMPUTER, "grfepc") ~ here::here("data", "2_results"),
  str_detect(COMPUTER, "epfl\\.ch") ~ here::here("data", "2_results"),
  TRUE ~ ""
)

DIR_FINEMAP <- glue::glue("{DIR_SCRATCH}/finemap")

## in case PATHOGEN are gene counts, load the length of the genes, so relative counts can be calculated
FILE_PATH_COUNTS <-
  glue::glue("{DIR_PATHOGEN}/NC_009334_with_t1_alts_gene_summary.max0samp.dat")

##////////////////////////////////////////////////////////////////
##                      define data                             //
##////////////////////////////////////////////////////////////////

## -- Define data to work with -------------------------------------------------

HOST <- "data/0_raw/host/geno/EBVG2G_05"

if (exists("PATHOGEN")) {
  ## extract name of pathogen (without path)
  NAM <- fs::path_file(PATHOGEN)
  
  ## same for GRM matrix
  PATHOGEN_FORGRM <- case_when(
    str_detect(PATHOGEN, "t1") ~ as.character(
      glue::glue(
        "{DIR_PATHOGEN}/NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat"
      )
    ),
    str_detect(PATHOGEN, "t2") ~ as.character(
      glue::glue(
        "{DIR_PATHOGEN}/NC_007605_with_t2_alts_aa_variant_matrix.non_synonymous.binary.dat"
      )
    )
  )
  
  
}


##////////////////////////////////////////////////////////////////
##                        colors                                //
##////////////////////////////////////////////////////////////////

cols <-
  colorblindr::palette_OkabeIto ## orange, bright blue, green, yellow, blue
# plot(1:8, col = colorblindr::palette_OkabeIto, pch = 16)

##////////////////////////////////////////////////////////////////
##                        debugging                             //
##////////////////////////////////////////////////////////////////


## set up host file
from_to_debugging <- 35642566 + c(-1, 1) * 200e3
chr_debugging <- 22


maf.thresh <- 0.05
hwe.thresh <- 1e-6
callrate.snps.thresh <- 0.01
callrate.inds.thresh <- 0.03
#callrate.snps.thresh <- 0.2 ## used to be 0.01, running into numerical problems with 0.2
#callrate.inds.thresh <- 0.4 ## used to be 0.03
# callrate.snps.thresh <- 0.1
# callrate.inds.thresh <- 0.1
# rueger@grfepc3:~/G2G-EBV/data/1_processed$ wc -l EBVG2G_SHCS_Sanger_QC.fam
# 268 EBVG2G_SHCS_Sanger_QC.fam
# rueger@grfepc3:~/G2G-EBV/data/1_processed$ wc -l EBVG2G_SHCS_Sanger_QC.bim
# 4431188 EBVG2G_SHCS_Sanger_QC.bim

## not used
## hz.thresh <- 3


##////////////////////////////////////////////////////////////////
##                          set bin                             //
##////////////////////////////////////////////////////////////////

PLINK <- here::here("bin", "plink2_linux")
PLINK1 <- here::here("bin", "plink1.9_linux")

## gcta:
## help: https://cnsgenomics.com/software/gcta/#ManipulatingtheGRM
GCTA <- here::here("bin", "gcta64_1.91.7_linux.beta")

## There is a option to provide multiple GRMs, but not sure if this only
## applies to chromosomes or really multiple random effects
#--mgrm multigrm.txt or --mgrm-bin multigrm.txt
#--mbfile chrs.txt


##////////////////////////////////////////////////////////////////
##                    extra functions                           //
##////////////////////////////////////////////////////////////////

## https://gist.github.com/DavisVaughan/f41027c5334f7debc09fc606798e20ef
str_ignore <- function(string, pattern) {
  string[!str_detect(string, pattern)]
}

p.THRESH.lenient <- 5e-8 ## checked with jacques


## sample size extraction from GCTA log files

sample_size_extract <-
  function(DIR = "/home/rueger/G2G-EBV/data/2_results", phenotype = "NC_009334_with_t1_alts_aa_variant_matrix.binary.dat")
  {
    files <-
      fs::dir_ls(DIR) %>% stringr::str_subset(phenotype) %>% stringr::str_subset("log")
    
    sample_size <-
      purrr::map(files, function(x)
        system(glue::glue("grep -F 'observations' {x}"), intern = TRUE) %>% stringr::str_split(" ") %>% purrr::map(1)) %>% unlist() %>% as.numeric()
    
    outcome_vec <-
      files %>% stringr::str_replace(DIR, "") %>% stringr::str_replace(glue::glue("/lmm_gcta_{phenotype}_"), "") %>% stringr::str_replace(".log", "")
    
    sample_size_df <-
      tibble::tibble(files = files,
                     n = sample_size,
                     outcome = outcome_vec)
    
    return(sample_size_df)
  }

#sample_size_extract(DIR = "/home/rueger/G2G-EBV/data/2_results", phenotype = "NC_009334_with_t1_alts_aa_variant_matrix.binary.dat")

## //////////////////
## PLOTs
## //////////////////
qqplot_wrapper <- function(x, title = "", subtitle = "")
{
  N <- length(x)
  
  ## calculate the expected axis
  expected <-
    sort(-log10((1:N) / N - 1 / (2 * N)))
  observed <-
    sort(-log10(x))
  
  qp <- ggplot() +
    geom_abline(intercept = 0,
                slope = 1,
                alpha = I(0.5)) +
    geom_point(aes(expected, observed), shape = 1) +
    #   geom_point(aes(expected, observed), shape = 16, alpha = I(0.3)) +
    xlab(expression(Expected ~ ~ -log[10](italic(p)))) +
    ylab(expression(Observed ~ ~ -log[10](italic(p)))) +
    labs(title = title, subtitle = subtitle)
  
  return(qp)
}

mhtplot_wrapper <-
  function(data , title, subtitle, THRESH1, THRESH2, color_highlight = "red") {
    ## chr
    ## pos
    ## p
    
    
    ## add small space
    ## equidistance
    data_max <-
      data %>% group_by(chr) %>% slice(which.max(pos)) %>% mutate(pos_max = TRUE)
    data <-
      left_join(data, data_max) %>% mutate(pos_max = case_when(is.na(pos_max) ~ FALSE, TRUE ~ pos_max))
    
    CONST <- 10000
    data <- data %>% dplyr::arrange(chr, pos) %>%
      dplyr::mutate(tmp = 1) %>%
      mutate(tmp = case_when(pos_max ~ tmp + CONST,
                             TRUE ~ tmp)) %>%
      mutate(cumsum.tmp = cumsum(tmp))
    
    ## real distance
    med.dat <-
      data %>% group_by(chr) %>% summarise(median.x = median(cumsum.tmp))
    
    ## signif snps
    CONST <- 4e5 ## same as in src/run_locuszoomplot
    data_snps_signif <-
      data %>% filter(p <= THRESH1) %>% select(chr, pos)
    data_snps_signif_range <-
      purrr::map2_dfr(data_snps_signif$chr, data_snps_signif$pos, function(x, y)
        data %>% filter(chr == x & pos < (y + CONST) & pos > (y - CONST)))
    
    
    if (nrow(data_snps_signif) > 0)
    {
      qp <- ggplot() +
        geom_point(data = data,
                   aes(cumsum.tmp,-log10(p), color = as.factor(chr)),
                   shape = 1) +
        #geom_point(data = data, aes(cumsum.tmp, -log10(p), color = as.factor(chr)),  shape = 16, alpha = I(0.5)) +
        labs(title = title) + xlab("Chromosomal Position") +
        ylab(expression(-log[10](italic(p)))) +
        scale_x_continuous(breaks = med.dat$median.x, labels = med.dat$chr) +
        geom_hline(yintercept = -log10(THRESH1),
                   linetype = 2) +
        geom_hline(
          yintercept = -log10(THRESH2),
          alpha = I(0.5),
          linetype = 2
        ) +
        scale_color_manual(values = rep(c("black", gray(0.4)), 11)) + theme(legend.position = "none") +
        geom_point(aes(cumsum.tmp,-log10(p)),
                   color = I(color_highlight),
                   data = data_snps_signif_range)
      
      
    } else{
      qp <- ggplot() +
        geom_point(data = data,
                   aes(cumsum.tmp,-log10(p), color = as.factor(chr)),
                   shape = 1) +
        #geom_point(data = data, aes(cumsum.tmp, -log10(p), color = as.factor(chr)),  shape = 16, alpha = I(0.5)) +
        labs(title = title) + xlab("Chromosomal Position") +
        ylab(expression(-log[10](italic(p)))) +
        scale_x_continuous(breaks = med.dat$median.x, labels = med.dat$chr) +
        geom_hline(yintercept = -log10(THRESH1),
                   linetype = 2) +
        geom_hline(
          yintercept = -log10(THRESH2),
          alpha = I(0.5),
          linetype = 2
        ) +
        scale_color_manual(values = rep(c("black", gray(0.4)), 11)) + theme(legend.position = "none")
      
    }
    return(qp)
    
    
  }









#' Find SNP identifiers based on Chromosome:Position
#'
#' @param chr numeric chromosome
#' @param pos b37
#' @param minor_allele the corresponding minor allele 
#' @param emsembl from snp.ensembl <- biomaRt::useEnsembl(biomart = "snp", dataset = "hsapiens_snp", GRCh = 37)
#'
#' @return character with refsnp_id

#'
#' @examples 
#' snp.ensembl <- biomaRt::useEnsembl(biomart = "snp", dataset = "hsapiens_snp", GRCh = 37)
#' chrpos2snp(chr = 8, pos = 35180595, ensembl = snp.ensembl) ## rs2950922
#' chrpos2snp(chr = 8, pos = 35180595, minor = "A", ensembl = snp.ensembl) ## rs2950922
#' chrpos2snp(chr = 8, pos = 35180595, minor = "G", ensembl = snp.ensembl) ## rs2950922


chrpos2snp <- function(chr, pos, minor_allele = NULL, ensembl){
  
  
  
  out <- biomaRt::getBM(
    attributes = c("refsnp_id", "allele"),#,  "chrom_start", "chrom_end" ,"allele"),
    filters = c("chr_name", "start", "end"), 
    values = list(chr, pos, pos), 
    mart = ensembl) 
  
  ## take alleles appart
  out <- out %>% separate(allele, c("A1", "A2"))
  
  
  if (is.null(minor_allele)) {
    return(out)
    
  }else {
    if(out$minor_allele == minor_allele){
      
      return(out %>% pull(refsnp_id))
      
    }else{
      warning(glue::glue("Minor alleles do not comply: ensemble minor allele is {out$minor_allele}, the data's minor allele is {minor_allele}")) 
    }
    
  }
  
}

