############################################################################
############################################################################
###                                                                      ###
###                           DESCRIPTION OF                             ###
###                            DATA USED                                 ###
###                                                                      ###
############################################################################
############################################################################

## this script calculates all the numbers needed to characterize the
## cohort data (SHCS and Sanger).

## make_4_descript.sh does the following:
## 1) create extra folder "data/4_descript"
## 2) run descript-make-gcta.R with a different DIR_PROCESSED for all 8 phenotype groups (Sanger and SHCS)
## 3) then read all the processed data in, plus all the raw one
## 4) stores output into DIR_DESCRIPT<- "/home/rueger/G2G-EBV/data/4_descript/"


## rmarkdown::render("src/descript-cohort.R")
## rstudioapi::viewer("src/descript-cohort.html")
## for technical stuff, see https://rmarkdown.rstudio.com/articles_report_from_r_script.html
## not working: knitr::opts_knit$set(results=FALSE, echo=FALSE)

## -------------------------------------------
##   measures of interest                   --
## -------------------------------------------

## - datasets
##     GT
##     PATHOGEN
##     COVAR
##     PC

## - QC procedure for each dataset
##     GT (MAF, etc)
##     PATHOGEN (frequency)
##     COVAR (-)
##     PC (-)

## - number of individuals pre- and post-QC
##     GT
##     PATHOGEN


## - number of genetic variants pre- and post-QC
##     GT
##     PATHOGEN

## - distribution of variables in each dataset
##     GT (SNPs histogram)
##     PATHOGEN (genetic variants, histogram)
##     COVAR (table)
##     PC (plot)


## ----------------------------------------------------------------
##   1. Setup                                                    --
## ----------------------------------------------------------------

## for which phenotype should summary be generated
## results_max3_bwa_mem_hiv_18_NC_009334_with_t1_alts_aa_variant_matrix.binary.dat
args <- c("/home/rueger/G2G-EBV/data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat")
   
args <- commandArgs(trailingOnly = TRUE)

PATHOGEN <- args[1]

## This R script creates all the variables
DIR_SRC <- here::here("src")
source(glue::glue("{DIR_SRC}/setup.R"))

## rename DIR_PROCESSED (this is where the data of this script lives)
NAM2 <-
  str_replace(PATHOGEN, "/home/rueger/G2G-EBV/data/0_raw/pathogen/", "") %>% str_replace("/", "_")
DIR_PROCESSED <-
  glue::glue("/home/rueger/G2G-EBV/data/4_descript/{NAM2}")
fs::dir_create(DIR_PROCESSED)

## Create short string for annotation purposes
NAM_short <- case_when(
  NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat" ~ "Rare EBV gene variation", ## binary
  NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_gene_matrix.max5samp.binary.non_synonymous.dat" ~ "Common EBV gene variation", ## binary
  NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat" ~ "EBV gene variation", ## continuous
  NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat" ~ "EBV amino acids" ## bianry)
)

## Create short string for annotation purposes
NAM_cohort <- case_when(
  str_detect(DIR_PROCESSED, "hiv") ~ "SHCS individuals",
  str_detect(DIR_PROCESSED, "sanger") ~ "GPC individuals"
)


type_outcome <- case_when(
  str_detect(NAM, "counts") ~ "continuous",
  TRUE ~ "binary"
)




## ----------------------------------------------------------------
##   2. Datasets                                                 --
## ----------------------------------------------------------------

## HOST ------------
HOST <-
  "/home/rueger/G2G-EBV/data/0_raw/host/geno/EBVG2G_SHCS_Sanger_rs"
## dpeends, see below at HOST_RAW

HOST_QC <-
  glue::glue("{HOST}_QC") %>% str_replace(DIR_HOST, DIR_PROCESSED)

fam_vec <- c("fid", "iid", "iid", "iid", "sex", "phenotype")
map_vec <- c("chr", "snpid", "pos_mor", "pos", "A1", "A2")

## this should be the very original dataset that Chris created before merging
if(str_detect(PATHOGEN, "hiv")) ## SHCS
{
  HOST_RAW <- glue::glue("{DIR_HOST}/EBVG2G_05")
}
if(str_detect(PATHOGEN, "sanger")) ## GPC
{
  HOST_RAW <- glue::glue("{DIR_HOST}/Sanger/merged")

}
fam_raw <- data.table::fread(glue::glue("{HOST_RAW}.fam"),
                    header = FALSE)
map_raw <- data.table::fread(glue::glue("{HOST_RAW}.bim"),
                             header = FALSE)

fam <- data.table::fread(glue::glue("{HOST_QC}.fam"),
                             header = FALSE)
map <- data.table::fread(glue::glue("{HOST_QC}.bim"),
                  header = FALSE)

names(fam) <- fam_vec
names(fam_raw) <- fam_vec
names(map) <- map_vec
names(map_raw) <- map_vec




## PATHOGEN ---------
phenotype_raw <- read_csv(PATHOGEN, na = "-9.0") %>%
  rename(id = sample) %>% janitor::clean_names()

phenotype <- read_delim(
  glue::glue(
    "{DIR_PROCESSED}/{dir(DIR_PROCESSED) %>% str_subset('pheno')}"
  ),
  col_names = TRUE,
  delim = " "
)


## COVAR ------------

covar_discrete <- read_delim(
  glue::glue("{DIR_PROCESSED}/forgcta_covar_discrete.covar"),
  col_names = TRUE,
  delim = " "
)

covar_numeric <- read_delim(
  glue::glue("{DIR_PROCESSED}/forgcta_covar_numeric.qcovar"),
  col_names = TRUE,
  delim = " "
)


plan <-
  read_delim(glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"), delim = " ")

## ----------------------------------------------------------------
##   3. QC                                                       --
## ----------------------------------------------------------------

## HOST
#system(
#  glue::glue(
#    "{PLINK} --bfile {HOST_QC}_intermediate --mind {callrate.inds.thresh} --maf {maf.thresh} --hwe {hwe.thresh} --make-bed --out {HOST_QC}"
#  )
#)

# plan <- plan %>% mutate(include = (outcome.freq >= outcome_thresh & outcome.freq <= (1-outcome_thresh)))


## ----------------------------------------------------------------
##   4. Number of individuals pre and post QC                    --
## ----------------------------------------------------------------
## - number of individuals pre-QC
##     GT
##     PATHOGEN

## - number of individuals post-QC
##     GT
##     PATHOGEN

## see below


## ----------------------------------------------------------------
##   5. Number of characteristics pre and post QC                --
## ----------------------------------------------------------------
## - number of genetic variants pre-QC
##     GT
##     PATHOGEN

## - number of genetic variants post-QC
##     GT
##     PATHOGEN

pathogen_text <- c(
  "Pathogen genome dataset:\n",
  glue::glue(
    "Before QC, {nrow(plan)} genetic variants were present and {nrow(phenotype_raw)} individuals.\n"
  ),
  glue::glue(
    "After QC (removing genetic variants with less than 10% frequency), {ncol(phenotype)-2} genetic variants were present and {nrow(phenotype)} individuals."
  ),
  "QC: outcome freq > 0.1"
)

host_text <- c(
  "Host genome variants dataset:\n",
  glue::glue(
    "Before QC, {nrow(map_raw)} genetic variants were present for {nrow(fam_raw)} individuals.\n"
  ),
  glue::glue(
    "After QC (see below), {nrow(map)} genetic variants were present for {nrow(fam)} individuals."
  ),
  glue::glue(
    "QC step 1, apply callrate SNPS: --geno {callrate.snps.thresh} --maf {maf.thresh} --hwe {hwe.thresh} --keep individuals-in-phenotype"
  ),
  glue::glue(
    "QC step 2, apply callrate individuals: --mind {callrate.inds.thresh} --maf {maf.thresh} --hwe {hwe.thresh}"
  )
)

text <- rbind(paste(pathogen_text, collapse = " "), paste(host_text, collapse = " ")) %>% as_tibble()

write_delim(text,
            path = glue::glue("{DIR_DESCRIPT}/descript_text_{NAM2}.txt"),
            delim = " ")

## ----------------------------------------------------------------
##   6. Feature distribution                                     --
## ----------------------------------------------------------------
## - distribution of variables in each dataset
##     GT (SNPs histogram)
##     PATHOGEN (genetic variants, histogram)
##     COVAR (table)
##     PC (plot)

## Distribution HOST --------------------------------------------

maf_host_raw <-
  gaston::read.bed.matrix(basename = HOST_RAW, rds = NULL)@snps
maf_host <-
  gaston::read.bed.matrix(basename = HOST_QC, rds = NULL)@snps
nrow_host <- format(nrow(maf_host), big.mark = "'")
nrow_host_raw <- format(nrow(maf_host_raw), big.mark = "'")

maf_host_all <-
  rbind(
    maf_host %>% select(id, maf) %>% mutate(type = glue::glue("post QC (n.variants = {nrow_host})")),
    maf_host_raw %>% select(id, maf) %>% mutate(type = glue::glue("pre QC (n.variants = {nrow_host_raw})"))
  ) %>% mutate(type = fct_relevel(factor(type), glue::glue("pre QC (n.variants = {nrow_host_raw})")))

qp_host <-
  ggplot(data = maf_host_all) + geom_histogram(aes(x = maf)) +
  xlab("MAF per variant") +
  labs(title = NAM_short,
       subtitle = "",
       caption = NAM_cohort) + facet_wrap( ~ type) + 
  xlim(0, 0.5)


## store
png(glue::glue("{DIR_DESCRIPT}/descript_host_{NAM2}.png"),
    height = 500, width = 1000)
print(qp_host)
dev.off()

## this is for figures-tables
write_delim(maf_host_all,
            path = glue::glue("{DIR_DESCRIPT}/descript_host_table_{NAM2}.txt"),
            delim = " ")


## Distribution PHENOTYPE --------------------------------------------
sm_pathogen <- function(x, type_outcome = c("binary", "continuous"))
{
  
  if (type_outcome == "binary") {
    out <- sum(x == 1, na.rm = TRUE) / length(x)
  }

  if (type_outcome == "continuous") {
    out <- mean(x, na.rm = TRUE)
  }
  return(out)
}


type_summary <- case_when(
  str_detect(NAM, "counts") ~ "Mean",
  TRUE ~ "Frequency"
)


maf_pathogen_raw <-
  apply(phenotype_raw %>% select(-1), 2, function(x)
    sm_pathogen(x, type_outcome = type_outcome)) %>% t() %>% as_tibble() %>% gather()

maf_pathogen <-
  apply(phenotype %>% select(-1,-2), 2, function(x)
    sm_pathogen(x, type_outcome = type_outcome)) %>% t() %>% as_tibble() %>% gather()

nrow_pathogen <- format(nrow(maf_pathogen), big.mark = "'")
nrow_pathogen_raw <- format(nrow(maf_pathogen_raw), big.mark = "'")

maf_pathogen_all <-
  rbind(
    maf_pathogen %>% mutate(type = glue::glue("post QC (n.variants = {nrow_pathogen})")),
    maf_pathogen_raw %>% mutate(type = glue::glue(
      "pre QC (n.variants = {nrow_pathogen_raw})"
    ))
  ) %>% 
  mutate(type = fct_relevel(factor(type), glue::glue("pre QC (n.variants = {nrow_pathogen_raw})")))

## this is for figures-tables.R
write_delim(maf_pathogen_all,
            path = glue::glue("{DIR_DESCRIPT}/descript_pathogen_table_{NAM2}.txt"),
            delim = " ")


qp_pathogen <-
  ggplot(data = maf_pathogen_all) + geom_histogram(aes(x = value)) +
  xlab(glue::glue("{type_summary} value per variant")) +
  labs(title = NAM_short,
       subtitle = "",
       caption = NAM_cohort) + facet_wrap( ~ type, scales = "free")


## store
png(glue::glue("{DIR_DESCRIPT}/descript_pathogen_{NAM2}.png"),
    height = 500, width = 1000)
print(qp_pathogen)
dev.off()



## Distribution PHENOTYPE per individual ---------------------------------------
maf_pathogen_per_id <-
  apply(phenotype %>% select(-1,-2), 1, function(x)
    sm_pathogen(x, type_outcome = type_outcome)) %>% t() %>% as_tibble() %>% gather()


## this is for figures-tables.R
write_delim(maf_pathogen_per_id,
            path = glue::glue("{DIR_DESCRIPT}/descript_pathogen_table_per_id_{NAM2}.txt"),
            delim = " ")


qp_pathogen_per_id <-
  ggplot(data = maf_pathogen_per_id) + geom_histogram(aes(x = value)) +
  xlab(glue::glue("{type_summary} value per variant")) +
  labs(title = NAM_short,
       subtitle = "",
       caption = NAM_cohort)

png(glue::glue("{DIR_DESCRIPT}/descript_pathogen_per_id_{NAM2}.png"),
    height = 500, width = 500)
print(qp_pathogen_per_id)
dev.off()


## Distribution covars --------------------------------------------

covars <- full_join(covar_numeric, covar_discrete)

if (all( c("age", "sex") %in% names(covars))){
  ## turn sex into a factor, and normalise all log10 variables again
  covars %<>% mutate(
    sex = forcats::fct_recode(as.factor(sex), female = "2", male = "1")
  ) %>% select(-X1,-X2)
  ##      Sex (1=male; 2=female; other=unknown)
  
  qp_covars <- GGally::ggpairs(covars %>% select(sex, age, type_1_2),
                               aes(colour = sex, alpha = 0.4)) #+ scale_color_manual(values=colorblindr::palette_OkabeIto[1:3])
  
  ## store
  png(
    glue::glue("{DIR_DESCRIPT}/descript_covars_{NAM2}.png"),
    width = 400,
    height = 400
  )
  print(qp_covars)
  dev.off()
  
}




## ----------------------------------------------------------------
##   6. Feature distribution TABLE                               --
## ----------------------------------------------------------------

## complete

## covars
if (all(covars %in% "sex" )){
  sm_covars <- covars %>% skimr::skim_to_wide() %>%
    rename(median = p50) %>%
    select(variable, n, complete, top_counts, mean, median, sd, p0, p100) %>%
    mutate(dataset = "Covariates") %>%
    rename(min = p0, max = p100)
 # mutate(top_counts = NA) %>%
    
}else{
  sm_covars <- covars %>% skimr::skim_to_wide() %>%
    rename(median = p50) %>%
    select(variable, n, complete, top_counts, mean, median, sd, p0, p100) %>%
    mutate(dataset = "Covariates") %>%
    rename(min = p0, max = p100)
  
}

## pathogen
sm_pathogen <- maf_pathogen_all %>%
  filter(str_detect(type, 'post QC')) %>%
  select(value) %>%
  # dplyr::rename(glue::glue(backtick("{type_summary} of present variant")) = value) %>%
  skimr::skim_to_wide() %>%
  rename(median = p50) %>%
  select(variable, mean, median, sd, p0, p100) %>%
  rename(min = p0, max = p100) %>%
  mutate(dataset = glue::glue("Pathogen genome ({nrow_pathogen} {NAM_short})"))

sm_pathogen$variable[1] <- glue::glue("Variant frequency")

if(type_outcome == "continuous") {
  sm_pathogen$variable[1] <- glue::glue("{type_summary} of variant fraction")
  
}

## host
sm_host <- maf_host_all %>%
  filter(str_detect(type, 'post QC')) %>%
  select(maf) %>% rename(`Minor allele frequency` = maf) %>%
  skimr::skim_to_wide() %>% rename(median = p50) %>%
  select(variable, mean, median, sd, p0, p100) %>%
  rename(min = p0, max = p100) %>%
  mutate(dataset = glue::glue("Host genome ({nrow_host} SNPs)"))

## rbind
out <-
  bind_rows(sm_host, sm_pathogen, sm_covars) %>% select(dataset,
                                                        variable,
                                                        n,
                                                        complete,
                                                        top_counts,
                                                        mean,
                                                        median,
                                                        sd,
                                                        min,
                                                        max) %>%
  mutate_at(.vars = vars(mean, median, sd, min, max), .funs = as.numeric) ## this is a bit risky, but does the job


write_delim(out,
            path = glue::glue("{DIR_DESCRIPT}/descript_table_{NAM2}.txt"),
            delim = " ")









## ----------------------------------------------------------------
##   7. Principal components                                     --
## ----------------------------------------------------------------





## ----------------------------------------------------------------
## HOST principal components                                 
## ----------------------------------------------------------------

## eigenvalues ---------
host_eigen <- read_delim(glue::glue("{DIR_PROCESSED}/host_grm_principal_components.eigenval"), delim = " ", col_names = FALSE)

qp_eigenvalue <- ggplot(data = host_eigen %>% mutate(x = 1:nrow(host_eigen))) + geom_point(aes(y = X1, x = x), shape = I(1)) + xlab("PC") + ylab("Eigenvalue")+ labs(title = "Eigenvalues of host principal components")
png(glue::glue("{DIR_DESCRIPT}/descript_host_pcs_eigenvalues_{NAM2}.png"), width = 500, height = 500)
qp_eigenvalue
dev.off()

## pcs -----------
host_pc <- read_delim(glue::glue("{DIR_PROCESSED}/host_grm_principal_components.eigenvec"), delim = " ", col_names = FALSE)
names(host_pc) <- c("X1", "X2", paste0("PC", 1:20))

p1 <- qplot(PC1, PC2, data = host_pc, shape = I(1)) + theme(legend.title = element_blank()) + labs(title = "Host principal components 1-4", subtitle = "PC-1 vs. PC-2")
p2 <- qplot(PC1, PC3, data = host_pc, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-1 vs. PC-3")
p3 <- qplot(PC1, PC4, data = host_pc, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-1 vs. PC-4")
p4 <- qplot(PC1, PC5, data = host_pc, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-1 vs. PC-5")
p5 <- qplot(PC1, PC6, data = host_pc, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-1 vs. PC-6")

p6 <- qplot(PC2, PC3, data = host_pc, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-2 vs. PC-3")
p7 <- qplot(PC2, PC4, data = host_pc, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-2 vs. PC-4")
p8 <- qplot(PC2, PC5, data = host_pc, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-2 vs. PC-5")
p9 <- qplot(PC2, PC6, data = host_pc, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-2 vs. PC-6")

p10 <- qplot(PC3, PC4, data = host_pc, shape = I(1) )+ labs(subtitle = "PC-3 vs. PC-4")
p11 <- qplot(PC3, PC5, data = host_pc, shape = I(1) )+ labs(subtitle = "PC-3 vs. PC-5")
p12 <- qplot(PC3, PC6, data = host_pc, shape = I(1) )+ labs(subtitle = "PC-3 vs. PC-6")

p13 <- qplot(PC4, PC5, data = host_pc, shape = I(1) )+ labs(subtitle = "PC-4 vs. PC-5")
p14 <- qplot(PC4, PC6, data = host_pc, shape = I(1) )+ labs(subtitle = "PC-4 vs. PC-6")

p15 <- qplot(PC5, PC6, data = host_pc, shape = I(1) )+ labs(caption = glue::glue("Displaying PCs of {nrow(host_pc)} individuals."), subtitle = "PC-5 vs. PC-6")


pempty <- qplot(data = NULL) + geom_blank() + theme_void()

library(patchwork)
qp <- p1 + p2 + p3 + p4 + p5 + pempty + p6 + p7 + p8 + p9 + pempty + pempty  + p10 +  p11 + p12+ pempty + pempty+ pempty + p13 + p14 + pempty + pempty +pempty + pempty + p15 + plot_layout(ncol = 5)#, tag_level = "new") + plot_annotation(tag_levels = c("A", "1"), tag_suffix = ")")
png(glue::glue("{DIR_DESCRIPT}/descript_host_pcs_{NAM2}.png"), width = 1200, height = 1340)
print(qp)
dev.off()


## ----------------------------------------------------------------
## PATHOGEN principal components                                 
## ----------------------------------------------------------------
## batch is apprently meaningless (I just have this memory, but I don't actually know)

## visualise batch
#PC <- PC_raw %>% select(PC1, PC2, PC3, PC4, id)# %>% left_join(dat_pathogen_batch %>% select(batch, X1), by = c("id" = "X1" ))
load(glue::glue("{DIR_PROCESSED}/forgcta_principal_components.RData"))

## eigenvalue -----------

PC_eigenvalues <- data.frame(PC_eigenvalues, x = 1:nrow(PC_eigenvalues))
qp_eigenvalue <- ggplot(data = PC_eigenvalues) + geom_point(aes(y = eigenvalue, x = x), shape = I(1)) + xlab("PC") + ylab("Eigenvalue")+ labs(title = "Eigenvalues of pathogen principal components")
png(glue::glue("{DIR_DESCRIPT}/descript_pathogen_pcs_eigenvalues_{NAM2}.png"), width = 500, height = 500)
qp_eigenvalue
dev.off()

## pcs -------------------

## remove individuals
PC <- PC %>% filter(id %in% host_pc$X1)


p1 <- qplot(PC1, PC2, data = PC, shape = I(1)) + theme(legend.title = element_blank())+ labs(title = "Pathogen principal components 1-4", subtitle = "PC-1 vs. PC-2")
p2 <- qplot(PC1, PC3, data = PC, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-1 vs. PC-3")
p3 <- qplot(PC1, PC4, data = PC, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-1 vs. PC-4")
p4 <- qplot(PC1, PC5, data = PC, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-1 vs. PC-5")
p5 <- qplot(PC1, PC6, data = PC, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-1 vs. PC-6")

p6 <- qplot(PC2, PC3, data = PC, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-2 vs. PC-3")
p7 <- qplot(PC2, PC4, data = PC, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-2 vs. PC-4")
p8 <- qplot(PC2, PC5, data = PC, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-2 vs. PC-5")
p9 <- qplot(PC2, PC6, data = PC, shape = I(1))+ theme(legend.title = element_blank())+ labs(subtitle = "PC-2 vs. PC-6")

p10 <- qplot(PC3, PC4, data = PC, shape = I(1) )+ labs(subtitle = "PC-3 vs. PC-4")
p11 <- qplot(PC3, PC5, data = PC, shape = I(1) )+ labs(subtitle = "PC-3 vs. PC-5")
p12 <- qplot(PC3, PC6, data = PC, shape = I(1) )+ labs(subtitle = "PC-3 vs. PC-6")

p13 <- qplot(PC4, PC5, data = PC, shape = I(1) )+ labs(subtitle = "PC-4 vs. PC-5")
p14 <- qplot(PC4, PC6, data = PC, shape = I(1) )+ labs(subtitle = "PC-4 vs. PC-6")

p15 <- qplot(PC5, PC6, data = PC, shape = I(1) )+ labs(caption = glue::glue("Displaying PCs of {nrow(PC)} individuals."), subtitle = "PC-5 vs. PC-6")


pempty <- qplot(data = NULL) + geom_blank() + theme_void()

library(patchwork)
qp <- p1 + p2 + p3 + p4 + p5 + pempty + p6 + p7 + p8 + p9 + pempty + pempty  + p10 +  p11 + p12+ pempty + pempty+ pempty + p13 + p14 + pempty + pempty +pempty + pempty + p15 + plot_layout(ncol = 5)#, tag_level = "new") + plot_annotation(tag_levels = c("A", "1"), tag_suffix = ")")

png(glue::glue("{DIR_DESCRIPT}/descript_pathogen_pcs_{NAM2}.png"), width = 1200, height = 1340)
print(qp)
dev.off()





