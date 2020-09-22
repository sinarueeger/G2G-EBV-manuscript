
###########################################################################
###########################################################################
###                                                                     ###
###                               RESULTS                               ###
###                                                                     ###
###########################################################################
###########################################################################

## Collect and visualize the results from data/2_results
##
## Make pretty plots

## 0. Setup
## 1. Quick & Dirty Search 
## 2. Genomic lambda
## 3. QQ plots and MHT plots


publication <- TRUE

##////////////////////////////////////////////////////////////////
##                           0. Setup                           //
##////////////////////////////////////////////////////////////////


#PATHOGEN=(/home/rueger/G2G-EBV/data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat ## 69
#          /home/rueger/G2G-EBV/data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat ##  572
## NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat

## some default args that will be overwritten in the next line, but here for debugging
args <-
  c(
    "data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat",
    "t1",
    "forreal",
    "gcta"
  )


args <- commandArgs(trailingOnly = TRUE)


## some default values
PATHOGEN <- args[1] ## 1:length(NAMS)
taxid <- args[2]
debugging <-
  ifelse(args[3] == "debugging", TRUE, FALSE) ## if "debugging" make datasets small for checking if it runs
random_outcome <-
  ifelse(args[3] == "random", TRUE, FALSE) ## if "random", create random outcome vectors with differing frequency of 1's and 0's
tool <- args[4] ## gcta or something else (like gaston)



## This R script creates all the variables
## takes counter_for_pathogen_outcome as input!
DIR_SRC <- here::here("src")
source(glue::glue("{DIR_SRC}/setup.R"))


## add string to NAM, to indicate debugging, random or validation mode
if (debugging)
  NAM <- glue::glue("{NAM}_debugging")
if (random_outcome)
  NAM <- glue::glue("{NAM}_random")

## default sample size
## CHECK in 2_results/*log
## run script in munge/samplesize.R for this
#N_GCTA <- 268 
N_GCTA_df <- sample_size_extract(DIR = DIR_SCRATCH, phenotype = NAM) ## looks for sample size in log files, stored in setup.R


## read in plan and add sample size
plan <-
  read_delim(glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"), delim = " ") %>% left_join(N_GCTA_df %>% select(outcome, n), by = "outcome")


## determine number of tests
n_tests <- 458 ## run script in src/eff-number-phenotypes.R for this



## determine significance threshold
if (debugging)
{
  THRESH <- 1
  
} else{
  
  THRESH <- p.THRESH.lenient / (n_tests) # 0.05/( 6365738 * n_tests)
  
}


## Create short string for annotation purposes
NAM_short <- case_when(
  NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat" ~ "Gene phenotype (binary, variants < 1 sample)",
  NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_gene_matrix.max5samp.binary.non_synonymous.dat" ~ "Gene phenotype (binary, variants < 5 samples )",
  NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat" ~ "Gene phenotype (variant frequency)",
  NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat" ~ "Amino acid phenotype"
)

NAM2 <- case_when(
  NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat" ~ "gene_matrix.max1samp",
  NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_gene_matrix.max5samp.binary.non_synonymous.dat" ~ "gene_matrix.max5samp",
  NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat" ~ "gene_matrix.max0samp",
  NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat" ~ "aa_variant_matrix"
)

##///////////////////////////////////////////////////////////////
##           Dig up set of SNPs with strict QC                 //
##///////////////////////////////////////////////////////////////

## done in make_2_results.sh (src/post_QC)

## Define HOST data
HOST <-
  plan %>% select(host) %>% distinct() %>% assertr::verify(nrow(.) == 1) %>% pull(host)
NAM_HOST <- fs::path_file(HOST)##  glue::glue("{(DIR_HOST)}/"))

## read the fam file (for later purposes)
strict_QC_snps <-
  read_delim(glue::glue("{DIR_SCRATCH}/{NAM_HOST}_QC.bim"),
             delim = "\t",
             col_names = FALSE) %>% pull(X2)



##///////////////////////////////////////////////////////////////
##                   1. Quick & Dirty Search                   //
##///////////////////////////////////////////////////////////////


## among all gwas summary stats, select for each GWAS only the snps that have P < 1e-8
system(glue::glue("./src/subset-results.sh {NAM} {DIR_SCRATCH} {debugging}")) ## this one is searching for P < 1e-8 (when no debugging), and P <1 when debugging
## kill all awk processes: kill $(pgrep -f wget)

if (debugging){
  ## read them all in
  dat <-
    fs::dir_ls(DIR_SCRATCH) %>%
    str_subset(NAM) %>%
    str_subset("subset.txt") %>%
    str_subset(tool) %>%
    str_subset("debugging") %>% 
    str_ignore("strictQC") %>%
    str_ignore("validation") %>%
    str_ignore("clumped") %>%
    str_ignore("log") %>%
    str_ignore("nosex") %>%
    purrr::map_dfr(
      function(x)
        fread(x, header = TRUE)  %>% mutate(
          global.outcome = NAM, 
          outcome = x%>% str_replace_all(
            glue::glue(
              "{DIR_SCRATCH}/|.mlma|.txt|{NAM}_|lmm_|{tool}_|-subset.txt"
            ),
            "") 
        )
    ) %>% 
    rename(Chr = V1, SNP = V2, bp = V3, effect.allele = V4, other.allele = V5, Freq = V6, b = V7, se = V8, p = V9)
  
} else {
  ## read them all in
  dat <-
    fs::dir_ls(DIR_SCRATCH) %>%
    str_subset(NAM) %>%
    str_subset("subset.txt") %>%
    str_subset(tool) %>%
    str_ignore("validation") %>%
    str_ignore("debugging") %>%
    str_ignore("strictQC") %>%
    str_ignore("clumped") %>%
    str_ignore("log") %>%
    str_ignore("nosex") %>%
    purrr::map_dfr(
      function(x)
        fread(x, header = FALSE)  %>% mutate(
          global.outcome = NAM, 
          outcome = x%>% str_replace_all(
            glue::glue(
              "{DIR_SCRATCH}/|.mlma|.txt|{NAM}_|lmm_|{tool}_|-subset.txt"
            ),
            "") 
        )
    ) %>% 
    rename(Chr = V1, SNP = V2, bp = V3, effect.allele = V4, other.allele = V5, Freq = V6, b = V7, se = V8, p = V9)
  
}
## dat %>% arrange(p) %>% slice(1:6)



## only consider good QC ones
dat <- dat %>% filter(SNP %in% strict_QC_snps)

## if all files are empty, stop the script
stopifnot(nrow(dat) > 0)


## select only subset that survives threshold
dat_select <-
  dat %>% filter(p <= THRESH) %>% arrange((p)) %>% left_join(plan %>% select(outcome, outcome.freq, n), by = "outcome")

dat_select_lenient <-
  dat %>% filter(p <= p.THRESH.lenient) %>% arrange((p)) %>% left_join(plan %>% select(outcome, outcome.freq, n), by = "outcome")


## if gene 1 or 5, only consider outcome.freq > 0.04
dat_select <- dat_select %>% filter(outcome.freq >= 0.04)
dat_select_lenient <- dat_select_lenient %>% filter(outcome.freq >= 0.04)


## add info about other gene type

if (FALSE)
{
  
  other.NAM <- case_when(
    NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat" ~ NA_character_,
    NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat" ~ "NC_009334_with_t1_alts_aa_variant_matrix.binary.dat",
    NAM %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_aa_variant_matrix.binary.dat" ~ "NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat"
  )
  
  dat_other <-
    purrr::map_dfr(unique(dat_select$outcome), function(out)
      fs::dir_ls(DIR_SCRATCH) %>%
        str_subset(other.NAM) %>%
        str_subset(out) %>%
        str_subset(tool) %>%
        str_ignore("subset") %>%
        purrr::map_dfr(function(x)
          fread(x)) %>%
        unique() %>%
        filter(id %in% dat_select$id))
  
  ## merge
  if (nrow(dat_other) > 0) {
    dat_select <-
      dat_select %>% left_join(
        dat_other %>% select(id, score, p, outcome, model, outcome_global),
        by = c("id", "model", "outcome"),
        suffix = c("", ".partner.signal")
      )
  }
}



##////////////////////////////////////////////////////////////////
##                  2. GENOMIC LAMBDA                           //
##////////////////////////////////////////////////////////////////

## define smallest set (unique outcome and global outcome).
## E.g. if there are multiple loci present within one outcome, 
## only one GC analysis or MHT plot is produced. 


sel <-
  dat_select %>% dplyr::select(global.outcome, outcome, outcome.freq) %>%
  unique() %>% mutate(lambda_gc = NA, sig_lambda_gc = FALSE)



## Loop through all possible results and compute genomic lambda

for (i in 1:nrow(sel)) {
  NAM <- sel[i, "global.outcome"]
  VARIANT <- sel[i, "outcome"]
  FREQ <- sel[i, "outcome.freq"]
  
  if (tool == "gcta")
    out <-
    data.table::fread(glue::glue("{DIR_SCRATCH}/lmm_gcta_{NAM}_{VARIANT}.mlma")) %>% rename(chr = Chr, pos = bp)
  
  ## remove bad QC snps
  out %<>% filter(SNP %in% strict_QC_snps)
  
  
  
  ## calculate GC lambda
  lambda_gc <- GWAS.utils::genomic_inflation(P = na.omit(out$p))
  sel$lambda_gc[i] <- lambda_gc
  
  ## still significant after correcting p-values for lambda_gc
  out <-
    out %>% mutate(
      Z = qchisq(p, 1, lower.tail = FALSE),
      Z_gc = Z / lambda_gc,
      p_gc = (1 - pchisq(Z_gc, 1))
    )
  
  sel$sig_lambda_gc[i] <- any(out$p_gc < THRESH)
}

# only if still remains
if (!debugging) {
  sel <- sel %>% filter(sig_lambda_gc)
}

## only write out what survives

dat_select <- dat_select %>% right_join(sel %>% select(-sig_lambda_gc))
dat_select_lenient <- dat_select_lenient %>% right_join(sel %>% select(-sig_lambda_gc))


## merge with rsids -------------------------

## rsids
rsids_ <- data.table::fread("/svshare/cthorball/Reference_files/ReplaceSnpID_1kg.list", header = FALSE)
rsids <- rsids_ %>% filter(V1 %in% out$SNP) 
chrpos2snp_using_file <- function(x) 
{
  ## add proper rsids
  out_rs <- x %>% 
    rename(SNP_chrpos = SNP) %>% 
    left_join(rsids, by = c("SNP_chrpos" = "V1")) %>% 
    rename(SNP = V2) %>% 
    mutate(SNP = case_when(
      !is.na(SNP) ~ SNP,
      is.na(SNP) ~ SNP_chrpos, 
      TRUE ~ NA_character_
    ))# %>%
#    select(-SNP_chrpos)
  return(out_rs)
  
}

dat_select <- chrpos2snp_using_file(dat_select)
dat_select_lenient <- chrpos2snp_using_file(dat_select_lenient)

## write out
write_delim(dat_select,
            path = glue::glue("{DIR_OUTPUT}/results_top_{NAM}.txt"))

write_delim(dat_select_lenient,
            path = glue::glue("{DIR_OUTPUT}/results_lenient_{NAM}.txt"))





##////////////////////////////////////////////////////////////////
##                         3a. Q-Q Plots                         //
##                         3b. MHT Plots                         //
##////////////////////////////////////////////////////////////////


## loop through all results to produce QQ and MHT plots

for (i in 1:nrow(sel)) {
  NAM <- sel[i, "global.outcome"]
  VARIANT <- sel[i, "outcome"]
  FREQ <- sel[i, "outcome.freq"]
  
  if (tool == "gcta")
    out <-
    data.table::fread(glue::glue("{DIR_SCRATCH}/lmm_gcta_{NAM}_{VARIANT}.mlma")) %>% rename(chr = Chr, pos = bp)
  
  ## remove bad QC snps
  out %<>% filter(SNP %in% strict_QC_snps)
  
  ## add proper rsids
  out <- chrpos2snp_using_file(out)

  ## store this out again
  data.table::fwrite(out, file = glue::glue("{DIR_OUTPUT}/lmm_gcta_{NAM}_{VARIANT}_strictQC.mlma"), sep = "\t")
  

  ##----------------------------
  ## store stuff for locuszoom
  ##----------------------------
  
  PATH_TO_FILE <- glue::glue("{DIR_OUTPUT}/lmm_gcta_{NAM}_{VARIANT}_strictQC_locuszoomplot.txt")
  
  out_locuszoom <- out %>% 
    mutate(SNP = case_when(
      str_detect(SNP, ":") ~ paste0("chr", chr, ":", pos),
      str_detect(SNP, "rs") ~ SNP,
      TRUE ~ SNP
      )
    ) %>% 
    select(SNP, p)
  data.table::fwrite(out_locuszoom, PATH_TO_FILE, sep = "\t") 
  
  
  TOP_SNP <-
    dat_select %>% filter(global.outcome == NAM & outcome == VARIANT) %>% 
    group_by(Chr) %>% 
    slice(which.min(p)) %>% 
    mutate(SNP = case_when(
      str_detect(SNP, ":") ~ paste0("chr", Chr, ":", bp),
      str_detect(SNP, "rs") ~ SNP,
      TRUE ~ SNP
    )
    ) %>%
    pull(SNP)
  
  for (k in 1:length(TOP_SNP))
    system(glue::glue("./src/run_locuszoom.sh {PATH_TO_FILE} {TOP_SNP[k]} {VARIANT}"))
  

  
  ##----------------------------
  ## Lambda QC adaptation
  ##----------------------------
  
  lambda_gc <- GWAS.utils::genomic_inflation(P = na.omit(out$p))# sel$lambda_gc[i]
  out <-
    out %>% mutate(
      Z = qchisq(p, 1, lower.tail = FALSE),
      Z_gc = Z / lambda_gc,
      p_gc = (1 - pchisq(Z_gc, 1))
    )
  
  
  ##----------------------------
  ##  Q-Q plot                --
  ##----------------------------
  
  
  # for publication
  if (publication){
    qp.qq <- qqplot_wrapper(x = out %>% filter(p_gc != 0) %>% pull(p), subtitle = glue::glue("GC lambda = {round(lambda_gc, 2)}"), title = glue::glue("{VARIANT}"))
    
  } else {
    # TODO: deal with p == 0
    #  %>% filter(p_gc != 0)
    
    qp.qq <- ggplot(out %>% filter(p_gc != 0)) +
      ggGWAS::stat_gwas_qq_hex(aes(y = p)) +
      geom_abline(intercept = 0, slope = 1) +
      labs(
        title = glue::glue("{NAM_short}"),
        subtitle = glue::glue("{VARIANT} (freq = {round(FREQ, 3)})"),
        caption = glue::glue(
          "{expression('fitted value of' ~ lambda)} = {round(lambda_gc, 2)}"
        )
      )
    
  }
  png(glue::glue(
    "{DIR_OUTPUT}/results_top_qqplot_{NAM}_{VARIANT}.png"
  ))
  print(qp.qq)
  dev.off()
  
  
  #ggsave(glue::glue("{DIR_OUTPUT}/results_top_qqplot_{NAM}_{VARIANT}.pdf"), qp.qq, device = cairo_pdf, height = 6, width = 6)
  
  ##----------------------------
  ##  MHT plot                --
  ##----------------------------
  
  png(glue::glue("{DIR_OUTPUT}/results_top_manhattan_{NAM}_{VARIANT}.png"),
      width = 900)
  if (debugging) {
    qqman::manhattan(out  %>% filter(p_gc != 0),
                     chr = "chr",
                     bp = "pos",
                     p = "p")
    
  } else {
    
    # for publication
    if (publication) {
      qp <-
        mhtplot_wrapper(
          data = out  %>% filter(p_gc != 0) %>% select(SNP, chr, pos, p),
          title = glue::glue("{VARIANT}"),
          THRESH1 = THRESH,
          THRESH2 = p.THRESH.lenient
        )
      print(qp)
    }else{
      qqman::manhattan(
        out  %>% filter(p_gc != 0 & p < 1e-4),
        chr = "chr",
        bp = "pos",
        p = "p"
      )
      
    }
    
    
  }
  dev.off()
  
  # if (!debugging)
  #  {
  #   qp <- mhtplot_wrapper(data = out  %>% filter(p_gc != 0) %>% select(SNP, chr, pos, p), title = glue::glue("{VARIANT} in SHCS"), THRESH1 = THRESH, THRESH2 = p.THRESH.lenient)
  #  ggsave(glue::glue("{DIR_OUTPUT}/results_top_manhattan_{NAM}_{VARIANT}.pdf"), qp, device = cairo_pdf, height = 6, width = 10)
  
  #}
  
  ##----------------------------
  ##  make locuszoom plot online
  ##----------------------------
  ## scp 
  
  ## center around rs7981812
  ##  lmm_gcta_NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat_BRLF1:p.Glu377Ala.mlma
  ##  lmm_gcta_NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat_validation_BRLF1:p.Glu377Ala.mlma
  ## also rs17017583 for GPC
  
  ## center around rs9815101
  ##  lmm_gcta_NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat_BCRF1.mlma
  ##  lmm_gcta_NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat_validation_BCRF1.mlma
  ## also rs11923671 for GPC
  
  ## stored in 3_output/locuszoomplots/
  
}



##////////////////////////////////////////////////////////////////
##                  4. FINEMAP INFO                             //
##////////////////////////////////////////////////////////////////

#z;ld;snp;config;cred;log;n_samples
#dataset1.z;dataset1.ld;dataset1.snp;dataset1.config;dataset1.cred;dataset1.log;5363
#dataset2.z;dataset2.ld;dataset2.snp;dataset2.config;dataset2.cred;dataset2.log;5363

window <- 2e6


sel_finemap <-
  dat_select %>% dplyr::select(n, Chr, bp, global.outcome, outcome, outcome.freq) %>%
  group_by(Chr, global.outcome, outcome, n) %>% 
  summarize(bp_min = min(bp), bp_max = max(bp), bp_median = median(bp)) %>% 
  mutate(bp_from = round(bp_median - window/2), bp_to = round(bp_median + window/2)) %>%
  mutate(bp_from = case_when(
    bp_from >= bp_min ~ as.double(bp_min),
    bp_from < bp_min ~ bp_from
  )) %>%
  mutate(bp_to = case_when(
    bp_to <= bp_max ~ as.double(bp_max),
    bp_to > bp_max ~ bp_to
  ))

if(NAM2 == "gene_matrix.max1samp" & any(dat_select$outcome == "BDLF1"))
{
  sel_finemap_additional <-
    dat_select_lenient %>% dplyr::select(n, Chr, bp, global.outcome, outcome, outcome.freq) %>%
    unique() %>%
    group_by(Chr, global.outcome, outcome, n) %>% 
    summarize(bp_min = min(bp), bp_max = max(bp), bp_median = median(bp)) %>% 
    mutate(bp_from = round(bp_median - window/2), bp_to = round(bp_median + window/2)) %>%
    mutate(bp_from = case_when(
      bp_from >= bp_min ~ as.double(bp_min),
      bp_from < bp_min ~ bp_from
    )) %>%
    mutate(bp_to = case_when(
      bp_to <= bp_max ~ as.double(bp_max),
      bp_to > bp_max ~ bp_to
    )) %>% 
    filter(outcome == "BDLF1")
  
  sel_finemap <- rbind(sel_finemap, sel_finemap_additional) %>% unique()
  
}


write_delim(
  sel_finemap %>% ungroup() %>% select(Chr),
  path = glue::glue("{DIR_FINEMAP}/loci_finemap_{NAM2}_chr"),
  delim = " ",
  col_names = FALSE
)
write_delim(
  sel_finemap %>% ungroup() %>% select(bp_from),
  path = glue::glue("{DIR_FINEMAP}/loci_finemap_{NAM2}_bp_from"),
  delim = " ",
  col_names = FALSE
)
write_delim(
  sel_finemap %>% ungroup() %>% select(bp_to),
  path = glue::glue("{DIR_FINEMAP}/loci_finemap_{NAM2}_bp_to"),
  delim = " ",
  col_names = FALSE
)

write_delim(
  sel_finemap %>% ungroup() %>% select(Chr, bp_from, bp_to, outcome, global.outcome),
  path = glue::glue("{DIR_FINEMAP}/loci_finemap_{NAM2}_fulldetails"),
  delim = " "
)


## so the bp_to and bp_from are the files from which we need to generate files

## create master

K <- 1:nrow(sel_finemap)
master <- data.frame(z = glue::glue("{DIR_FINEMAP}/locus-{NAM2}-{K}.z"),
                     ld =glue::glue("{DIR_FINEMAP}/locus-{NAM2}-{K}.ld"),
                     snp =glue::glue("{DIR_FINEMAP}/locus-{NAM2}-{K}.snp"),
                     config =glue::glue("{DIR_FINEMAP}/locus-{NAM2}-{K}.config"),
                     cred =glue::glue("{DIR_FINEMAP}/locus-{NAM2}-{K}.cred"),
                     log =glue::glue("{DIR_FINEMAP}/locus-{NAM2}-{K}.log"),
                     n_samples = sel_finemap$n
)

write_delim(master, path = glue::glue("{DIR_FINEMAP}/master_{NAM2}"), delim = ";")

## create locus-$k.z files and locus incl variants
bim <- data.table::fread(glue::glue("{DIR_FINEMAP}/EBVG2G_05_QC_ldstore.bim"), header = FALSE)

for (i in 1:nrow(sel_finemap)) {
  NAM <- sel_finemap[i, "global.outcome"]
  VARIANT <- sel_finemap[i, "outcome"]
  
  out <-
    data.table::fread(glue::glue("{DIR_SCRATCH}/lmm_gcta_{NAM}_{VARIANT}.mlma")) %>% 
    rename(chr = Chr, pos = bp)
  
  ## remove bad QC snps
  out %<>% filter(SNP %in% strict_QC_snps)
  
  
  ## 
  data_z <- out %>% 
    filter(chr == sel_finemap$Chr[i] & dplyr::between(pos, sel_finemap$bp_from[i], sel_finemap$bp_to[i]) ) %>%
    select(SNP, chr, pos, A1, A2, Freq, b, se) %>%
    rename(rsid = SNP, chromosome = chr, position = pos, allele1 = A1, allele2 = A2, maf = Freq, beta = b, se = se) %>%
    mutate(maf = GWAS.utils::eaf2maf(maf))
  
  ## check if they are in the bim file
  data_z <- data_z %>% filter(rsid %in% bim$V2)
  
  write_delim(data_z, 
              path = glue::glue("{DIR_FINEMAP}/locus-{NAM2}-{i}.z"), 
              delim = " ")
  
  ## B)
  # create data.z
  #rsid chromosome position allele1 allele2 maf beta se
  #rs1         10 1 T C 0.35 0.0050 0.0208
  #rs2         10 1 A G 0.04 0.0368 0.0761
  #rs3         10 1 G A 0.18 0.0228 0.0199
  
  #SNP, Chr, bp, A1, A2, Freq, b, se
  #Chr     SNP     bp      A1      A2      Freq    b       se      p
  #1       1:11008:C:G     11008   G       C       0.0769231       4.98124e-05     0.000256642     0.846103
  #1       1:11012:C:G     11012   G       C       0.0769231       4.98124e-05     0.000256642     0.846103
  
  ## C)
  ## create ${DIR_PROCESSED}/locus-$LOCUS-incl-variants
  ## --incl-variants 		Extract LD information for variants given in the specified text file. 
  ## The specified file has 5 columns with a header: RSID, position, chromosome, A_allele and B_allele 		Requires --matrix or --table
  incl_variants <- data_z %>% 
    select(rsid, position, chromosome, allele1, allele2) %>%
    rename(RSID = rsid, A_allele = allele1, B_allele = allele2)
  
  write_delim(incl_variants, 
              path = glue::glue("{DIR_FINEMAP}/locus-{NAM2}-{i}-incl-variants"), 
              delim = " ")
  
}

system(glue::glue("./src/run_finemap.sh {NAM2}"))


## results of finemap
## --------------------------

#dat_select <- read_delim(glue::glue("{DIR_OUTPUT}/results_top_{NAM}.txt"), delim = " ")
#dat_select_lenient <- read_delim(glue::glue("{DIR_OUTPUT}/results_lenient_{NAM}.txt"), delim = " ")


dat_select$causal_finemap <- FALSE
dat_select$causal_finemap_prob <- NA
dat_select_lenient$causal_finemap <- FALSE
dat_select_lenient$causal_finemap_prob <- NA
res_tmp <- NULL

for (LOCUS in 1:nrow(sel_finemap))
{
  ##loop over all K's
  config <-
    data.table::fread(glue::glue("data/2_results/finemap/locus-{NAM2}-{LOCUS}.config"))
  #snp <- data.table::fread(glue::glue("data/2_results/finemap/locus-{NAM2}-{LOCUS}.snp"))
  ## >> the important bit (I believe) credibility set
  #cred <- data.table::fread(glue::glue("data/2_results/finemap/locus-{NAM2}-{LOCUS}.cred"))
  
  ## add to sel_finemap_results
  res_tmp <- rbind(res_tmp, config[1,]) 
  
  ## add to dat_select & dat_select_lenient
  causal_variant <-
    config %>% filter(rank == 1) %>% select(config) %>% pull(config) %>% str_split(",") %>% unlist()
  
  prob <- config %>% filter(rank == 1) %>% select(prob) %>% pull(prob)

  ## sets the causal vairant to TRUE, all other variants are FALSE
  dat_select$causal_finemap[dat_select$SNP_chrpos %in% causal_variant &
                              dat_select$outcome == sel_finemap$outcome[LOCUS]] <-
    TRUE
  dat_select$causal_finemap_prob[dat_select$SNP_chrpos %in% causal_variant &
                                   dat_select$outcome == sel_finemap$outcome[LOCUS]] <-
    prob
  
  dat_select_lenient$causal_finemap[dat_select_lenient$SNP_chrpos %in% causal_variant &
                                      dat_select_lenient$outcome == sel_finemap$outcome[LOCUS]] <-
    TRUE
  dat_select_lenient$causal_finemap_prob[dat_select_lenient$SNP_chrpos %in% causal_variant &
                                           dat_select_lenient$outcome == sel_finemap$outcome[LOCUS]] <-
    prob
  
}

## write out
sel_finemap_results <- cbind(sel_finemap %>% ungroup(), res_tmp %>% distinct())

write_delim(sel_finemap_results,
            path = glue::glue("{DIR_OUTPUT}/results_finemap_{NAM}.txt"))

write_delim(dat_select,
            path = glue::glue("{DIR_OUTPUT}/results_top_{NAM}.txt"))

write_delim(dat_select_lenient,
            path = glue::glue("{DIR_OUTPUT}/results_lenient_{NAM}.txt"))


##////////////////////////////////////////////////////////////////
##                  5. VIRALLOAD GWAS                           //
##////////////////////////////////////////////////////////////////

out <- data.table::fread(glue::glue("{DIR_SCRATCH}/GWAS_viralload.mlma"))

out %>% arrange(p) %>% head()

out %>% filter(p < 5e-8)

ggplot(out ) +
  ggGWAS::stat_gwas_qq_hex(aes(y = p)) +
  geom_abline(intercept = 0, slope = 1)



##////////////////////////////////////////////////////////////////
##               6. VIRALLOAD pathogen GWAS                     //
##////////////////////////////////////////////////////////////////

## This R script creates all the variables
DIR_SRC <- here::here("src")
source(glue::glue("{DIR_SRC}/setup.R"))

n_tests <- 458 ## run script in src/eff-number-phenotypes.R for this

pgwas <-
  fs::dir_ls(DIR_SCRATCH) %>%
  str_subset("pGWAS") %>%
  purrr::map_dfr(
    function(x)
      read_csv(x)  %>% mutate(
        outcome = x %>% str_replace_all(
          glue::glue(
            "{DIR_SCRATCH}/|.csv"
          ),
          "") 
      )
  ) %>% arrange(p.value)

## 767 tests, but nothing passes the threshold of 0.05/n_tests = 0.0001091703. 
pgwas %>% 
  filter(p.value < (0.05/n_tests))

ggplot(pgwas) + ggGWAS::stat_gwas_qq(aes(y = p.value, group = outcome, color = outcome)) + 
  geom_abline(intercept = 0, slope = 1)

