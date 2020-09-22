###########################################################################
###########################################################################
###                                                                     ###
###                         FIGURES + TABLES                            ###
###                                                                     ###
###########################################################################
###########################################################################
## making figures for manuscript
###########################################################################



## what we want to produce

## Table 1: Summary of Covariats, pathogen variants and human variants
## Figure 1: Descriptive Covariates
## Figure 2: MAF distribution of pathogen variants 
## Figure 3: MAF distribution of human variants
## >>> all in report/descript-cohort.R

## Figure 4: PCs (heatmap, batch)
## Table 2: results


## -------------------------------------------
##   set paths                              --
## -------------------------------------------

## >>> CHECK all results in 4_descript
PATHOGEN_GENE_BINARY_1 <- "/home/rueger/G2G-EBV/data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat"
PATHOGEN_GENE_BINARY_5 <- "/home/rueger/G2G-EBV/data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max5samp.binary.non_synonymous.dat"
PATHOGEN_GENE <- "/home/rueger/G2G-EBV/data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat"
PATHOGEN_AA <- "/home/rueger/G2G-EBV/data/0_raw/pathogen/results_max3_bwa_mem_hiv_18/NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat"
## they both have for host data etc the same values, checked in 4_descript

## This R script creates all the variables
DIR_SRC <- here::here("src")
source(glue::glue("{DIR_SRC}/setup.R"))


## used to read from the correct file
NAM_GENE_BINARY_1 <-
  str_replace(PATHOGEN_GENE_BINARY_1, "/home/rueger/G2G-EBV/data/0_raw/pathogen/", "") %>% str_replace("/", "_") 
NAM_GENE_BINARY_5 <-
  str_replace(PATHOGEN_GENE_BINARY_5, "/home/rueger/G2G-EBV/data/0_raw/pathogen/", "") %>% str_replace("/", "_") 
NAM_GENE <-
  str_replace(PATHOGEN_GENE, "/home/rueger/G2G-EBV/data/0_raw/pathogen/", "") %>% str_replace("/", "_") 
NAM_AA <-
  str_replace(PATHOGEN_AA, "/home/rueger/G2G-EBV/data/0_raw/pathogen/", "") %>% str_replace("/", "_") 

## used to read from the correct folder
DIR_PROCESSED_SHCS_GENE <- glue::glue("/home/rueger/G2G-EBV/data/4_descript/{NAM_GENE_BINARY_1}")
stopifnot(fs::dir_exists(DIR_PROCESSED_SHCS_GENE))


## -------------------------------------------
##  Plot Phenotypes                         --
## -------------------------------------------
maf_ebv <-
  fs::dir_ls(DIR_DESCRIPT) %>%
  str_subset("txt") %>% 
  str_ignore("sanger") %>% 
  str_ignore("per_id") %>% 
  str_subset("pathogen") %>% 
  str_subset("max1samp|aa_variant") %>%
  purrr::map_dfr(function(x)
    readr::read_delim(x, delim = " ") %>%
      filter(str_detect(type, "post")) %>%
      mutate(type = x)) %>% mutate(
        type = as.factor(type)) %>% mutate(type = fct_recode(type, 
                                                             `EBV amino acids` = "/home/rueger/G2G-EBV/data/4_descript/descript_pathogen_table_results_max3_bwa_mem_hiv_18_NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat.txt",
                                                             `EBV genes` = "/home/rueger/G2G-EBV/data/4_descript/descript_pathogen_table_results_max3_bwa_mem_hiv_18_NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat.txt",
        )
      )


## remove the ones with less than 0.04
maf_ebv <- maf_ebv %>% filter(value >= 0.04)

## add how many points
maf_ebv <- right_join(maf_ebv, maf_ebv %>% group_by(type) %>% summarize(n=n()) %>% mutate(dataset = glue::glue("{type} (#={n})")))


qp_ebv <- ggplot(data = maf_ebv, aes(value, group = type)) + geom_histogram(aes(y=..density..), color = gray(0.5), alpha=0.5, 
                                                position="identity") + xlim(c(0, NA)) +  
  guides(fill=FALSE, color = FALSE) + 
  facet_wrap(~dataset, scale = "free") + ## adding scale free x did not work
  labs(title = "Frequency distribution of EBV outcomes grouped by dataset",
     subtitle = glue::glue("{nrow(maf_ebv)} outcomes"),
     caption = "", 
     x = "Frequency") 

png(glue::glue("{DIR_DESCRIPT}/descript_pathogen.png"), height = 350, width = 600)
print(qp_ebv)
dev.off()


## -------------------------------------------
##  Plot Phenotypes (by ID)                 --
## -------------------------------------------
maf_ebv_id <-
  fs::dir_ls(DIR_DESCRIPT) %>% str_subset("txt") %>% 
  str_ignore("sanger") %>% 
  str_ignore("pathogen_table_results") %>% 
  str_subset("pathogen") %>% 
  str_subset("max1samp|aa_variant") %>%
  purrr::map_dfr(function(x)
    readr::read_delim(x, delim = " ") %>%
      mutate(type = x)) %>% mutate(
        type = as.factor(type)) %>% mutate(type = fct_recode(type, 
                                                             `EBV amino acids` = "/home/rueger/G2G-EBV/data/4_descript/descript_pathogen_table_per_id_results_max3_bwa_mem_hiv_18_NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat.txt",
                                                             `EBV genes` = "/home/rueger/G2G-EBV/data/4_descript/descript_pathogen_table_per_id_results_max3_bwa_mem_hiv_18_NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat.txt",
        )
        )

## remove the ones with less than 0.04
maf_ebv_id <- maf_ebv_id %>% filter(value >= 0.04)


## add how many points
maf_ebv_id <- right_join(maf_ebv_id, maf_ebv_id %>% group_by(type) %>% summarize(n=n()) %>% mutate(dataset = glue::glue("{type} (#={n})")))

qp_ebv_id <- ggplot(data = maf_ebv_id, aes(value, group = type)) + geom_histogram(aes(y=..density..), color = gray(0.5), alpha=0.5, 
                                                                                             position="identity") + xlim(c(0, NA)) + 
   guides(fill=FALSE, color = FALSE) + 
  facet_wrap(~dataset, scale = "free") + ## adding scale free x did not work
  labs(title = "EBV outcome distribution per individuals grouped by dataset",
       subtitle = glue::glue("{nrow(maf_ebv_id)} outcomes"),
       caption = "", 
       x = "Frequency") 




png(glue::glue("{DIR_DESCRIPT}/descript_pathogen_id.png"), height = 350, width = 600)
print(qp_ebv_id)
dev.off()


## -------------------------------------------
##  Plot Hostdata                           --
## -------------------------------------------

maf_host <- readr::read_delim(glue::glue("{DIR_DESCRIPT}/descript_host_table_{NAM_GENE_BINARY_1}.txt"), delim = " ") %>% filter(type == "post QC (n.variants = 4'291'179)")

qp_host_qc <-
  ggplot(data = maf_host) + geom_histogram(aes(x = maf, y=..density..), alpha = 0.8, position="identity") + geom_density(aes(x = maf)) + 
  xlab("Minor allele frequency") +
  ylab("Frequency") + 
  xlim(c(0, 0.5)) +
  labs(title = "Frequency distribution of human SNPs",
       caption = "",
       subtitle = glue::glue("After QC (# SNPs = {nrow(maf_host)})")) #+ facet_wrap( ~ type)

## store
png(glue::glue("{DIR_DESCRIPT}/descript_host.png"),
    height = 400, width = 500)
print(qp_host_qc)
dev.off()



## -------------------------------------------
##  Table results
## -------------------------------------------



## determine number of tests
n_tests <- 458 ## run script in munge/eff-number-phenotypes.R for this
THRESH <- p.THRESH.lenient / (n_tests) # 0.05/( 6365738 * n_tests)

## merging the two tables (SHCS and GPC) -----------------------------------

## SHCS
## b = log(or) >> http://cnsgenomics.com/software/gcta/GCTA_UserManual_v1.24.pdf
res <- fs::dir_ls(DIR_OUTPUT) %>%
  str_subset("results_top_") %>%
  str_ignore("validation") %>%
  str_subset(".txt") %>%
  purrr::map_dfr(function(x)
    fread(x)) %>%
  mutate(OR = case_when(
                 global.outcome == "NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat" ~ NA_real_,
                 TRUE ~ exp(b)
               )) %>% 
  dplyr::filter(p < THRESH) 

## replace global.outcome
res <- res %>% mutate(global.outcome = case_when(
  str_detect(global.outcome, "gene_matrix.max0samp.depth") ~ "gene (frequency)",
  str_detect(global.outcome, "_aa_") ~ "amino acid (binary)",
  str_detect(global.outcome, "gene_matrix.max1samp.binary") ~ "gene (binary, variants < 1 sample)",
  str_detect(global.outcome, "gene_matrix.max5samp.binary") ~ "gene (binary, variants < 5 samples)"
))


## define clumps --------------------------

## sanity check (make sure all loci not urther than {range} away)
range <- 1e6 ## if SNPs further away from that, define separate clumps
tmp <- res %>% group_by(Chr, outcome) %>% mutate(bp.diff = abs(bp - dplyr::lag(bp))) %>% dplyr::select(Chr, outcome, bp.diff)
stopifnot(all(tmp$bp.diff < range, na.rm = TRUE))

clumps <- res %>% dplyr::select(Chr, outcome) %>% distinct() %>%
  mutate(loci = 1:n())

res <- res %>% 
  left_join(clumps) %>% 
  dplyr::select(loci, everything())


## define top snp -------------------------


res_min <- res %>% group_by(loci) %>% slice(which.min(p)) %>% dplyr::select(loci, outcome, Chr, SNP) %>% mutate(top_SNP = TRUE)
res <- res %>% left_join(res_min) %>% mutate(top_SNP = replace_na(top_SNP, FALSE))


## add temporary SNP id for chr:pos ---------

## check if really no SNPid (the eugene eqtl function only works with SNPids)
res_no_rs <- res %>% filter(str_detect(SNP, ":"))


if (nrow(res_no_rs) > 0)
{
  snp.ensembl <-
    biomaRt::useEnsembl(biomart = "snp",
                        dataset = "hsapiens_snp",
                        GRCh = 37)
  
  plan(multiprocess)
  snps <-
    furrr::future_map2_dfr(res_no_rs$Chr, res_no_rs$bp, function(.x, .y) {
      out <-
        chrpos2snp(chr = .x,
                   pos = .y,
                   ensembl = snp.ensembl)
      
      
      if (nrow(out) > 0) {
        return(cbind(out, Chr = .x, pos = .y))
        
      }
      
    })
  
  ## merge SNPs
  ## doing this separately, bc joining by alleles won't work.
  ## and for del joining by alleles is crucial.
  res_snp <-
    left_join(
      res %>% filter(nchar(effect.allele) <= 1 &
                       nchar(other.allele) <= 1),
      snps %>% rename(SNP_tmp = refsnp_id, bp = pos) %>% select(SNP_tmp, Chr, bp, A1, A2),
      by = c("Chr", "bp", "other.allele" = "A1")
    )
  
  ## merge deletions (kinda get it to work by looking at the res and snps object)
  res_indel <-
    left_join(
      res %>% filter(nchar(effect.allele) > 1 |
                       nchar(other.allele) > 1),
      snps %>% rename(SNP_tmp = refsnp_id, bp = pos) %>% select(SNP_tmp, Chr, bp, A1, A2),
      by = c("Chr", "bp", "effect.allele" = "A1")
    )
  
  ## get it back together
  res <- rbind(res_snp, res_indel) %>% select(-A2)
  
  res <- res %>% mutate(SNP_tmp = case_when(is.na(SNP_tmp) ~ SNP,
                                            TRUE ~ SNP_tmp))
  
}

## LD with other SNPs --------------------------
# 
# f.extract.ld <- function(SNP.id2 = NULL, SNP.id1 = NULL, POP = NULL)
# {
#   # ext <- glue::glue("/ld/human/pairwise/{SNP.id1}/{SNP.id2}/1000GENOMES:phase_3:{POP}")  ## 
#   ext <- glue::glue("/ld/human/pairwise/{SNP.id1}/{SNP.id2}/") 
#   
#   server <- "https://rest.ensembl.org"
#   
#   r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
#   stop_for_status(r)
#   
#   out <- as_tibble(fromJSON(toJSON(content(r)))) %>% unnest() %>% filter(stringr::str_detect(population_name, POP))
#   return(out)
# }
# 
# 
# for (i in 5:7)#unique(res$loci))
# {
#   snp <- res %>% filter(loci == i & top_SNP) %>% pull(SNP)
#   other.snps <- res %>% filter(loci == i & !top_SNP) %>% filter(str_detect(SNP, "rs")) %>% pull(SNP)
#   
#   if (length(other.snps) > 0)
#   {
#     LD.SNP.SNPs <- purrr::map_df(other.snps, f.extract.ld, snp, "EUR") %>% mutate(r2 = as.numeric(r2)) %>% bind_rows() %>% unnest()
#     print(LD.SNP.SNPs)
#     
#   }  
#   
# }






## add gene, function, eqtl
## ----------------------------

#- eQTL
#- gene
#- gene description EBV

#library(biomaRt)
gene.ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37) # we will need an additional mart for genes
snp.ensembl <- biomaRt::useEnsembl(biomart = "snp", dataset = "hsapiens_snp", GRCh = 37)
biomaRt::listAttributes(snp.ensembl)
biomaRt::listAttributes(gene.ensembl)

cis_raw <- read_delim("/home/rueger/EQTL_DB/cis-eQTL_significant_20181017.txt.gz", delim = "\t") 
trans_raw <- read_delim("/home/rueger/EQTL_DB/trans-eQTL_significant_20181017.txt.gz", delim = "\t") 

source(glue::glue("{DIR_SRC}/function-eqtl.R"))


SNP <- res$SNP_tmp[1]
annotation_raw_gene <- parallel::mclapply(1:nrow(res), function(k)
{
  
  cat("START", k, "\n")
  SNP <- res$SNP[k]
  chr <- res$Chr[k]
  pos <- res$bp[k]
  
  
  ## only look for SNPs
  out_bm <- data.frame()
  
  if (str_detect(SNP, "rs"))
  {
    out_bm <- biomaRt::getBM(
      attributes = c(
        "ensembl_gene_stable_id",
        "refsnp_id",
        "chr_name",
        #       "ensembl_transcript_stable_id",
        "consequence_type_tv",
        "clinical_significance",
        "polyphen_score",
        "sift_score",
        "associated_gene",
        "phenotype_description",
        "variation_names"
      ),
      filters = "snp_filter", 
      values = SNP, 
      mart = snp.ensembl
    ) %>% slice(1) ####!!! careful here
    
  # else { ## if chr pos
    # out_bm <- biomaRt::getBM(
    #   attributes = c(
    #     "ensembl_gene_stable_id",
    #     "refsnp_id",
    #     "chr_name",
    #     #       "ensembl_transcript_stable_id",
    #     "consequence_type_tv",
    #     "clinical_significance",
    #     "polyphen_score",
    #     "sift_score",
    #     "associated_gene",
    #     "phenotype_description",
    #     "variation_names"
    #   ),
    #   filters = c("chr_name", "start", "end"), 
    #   values = list(chr, pos, pos), 
    #   mart = snp.ensembl
    # ) %>% slice(1) ####!!! careful here
    
  #}
    
    ## adding GENE -----------------------------------
    ## because we are using positions from GRCh = 37 in a next query, we need to pass that information on.
    out_bm_gene <-
      biomaRt::getBM(
        attributes = c('external_gene_name', 'ensembl_gene_id'),
        filters = c('ensembl_gene_id'),
        values = unique(out_bm$ensembl_gene_stable_id),
        mart = gene.ensembl
      )
    
    if (nrow(out_bm_gene) > 0)
    {
      out_bm <-
        out_bm %>% full_join(out_bm_gene,
                             by = c("ensembl_gene_stable_id" = "ensembl_gene_id"))
    }
  }  
  
  return(out_bm)
}, mc.cores = 12)
## won't work: 1:160072420:T:TCCAC >> in gene AL121987.1
## won't work: 21:38200666:GGAATGAAT:G


## gtex -----------------------------------
annotation_raw_gtex <-
  parallel::mclapply(1:nrow(res), function(k)
  {
    cat("START", k, "\n")
    SNP <- res$SNP[k]
    chr <- res$Chr[k]
    pos <- res$bp[k]
    
    out_gtex_files <-
      extract_gtex_files(chr = chr, pos = pos, rsid = SNP)
    
    return(out_gtex_files)
  }, mc.cores = 12)


## eugene ------------------------------------
annotation_raw_eugene <-
  parallel::mclapply(1:nrow(res), function(k)
  {
    cat("START", k, "\n")
    SNP <- res$SNP[k]
    chr <- res$Chr[k]
    pos <- res$bp[k]
    
    out_eugene <- extract_eugene(rsid = SNP)
    
    return(out_eugene)
  }, mc.cores = 12)

## eqtlgen ------------------------------------
annotation_raw_eqtlgen <-
  parallel::mclapply(1:nrow(res), function(k)
  {
    cat("START", k, "\n")
    SNP <- res$SNP[k]
    chr <- res$Chr[k]
    pos <- res$bp[k]
    
    out_eqtlgen <-
      extract_eqtlgen(
        chr = chr,
        pos = pos,
        rsid = SNP,
        cis_summarystats = cis_raw,
        trans_summarystats = trans_raw
      )
    
    
    return(out_eqtlgen)
  }, mc.cores = 12)


    

## form annotation data frame ------------------------------
## only include eqtls
#stopifnot(ncol(annotation_raw %>% bind_rows) == 13) ## length(annotation_raw %>% bind_rows) ## if more columns/anything changes, change the code below!

annotation_gene <- annotation_raw_gene %>% 
  bind_rows() %>% 
  discard(~all(is.na(.x)))
  
annotation_eqtlgen <- annotation_raw_eqtlgen %>% 
  bind_rows() %>% 
  group_by(SNP) %>%
  arrange(cis_eqtlgen.pvalue) %>% 
  summarize(cis_eqtlgen = paste(glue::glue("{cis_eqtlgen.gene} in {eqtlgen.tissue} (P={cis_eqtlgen.pvalue})"), collapse = ", ") ) %>%
#  mutate(cis_eqtlgen = case_when(cis_eqtlgen == "NA (P=NA)" ~ NA_character_,
#                              TRUE ~ cis_eqtlgen)) %>% 
  ungroup()


annotation_eugene <- annotation_raw_eugene %>% 
  bind_rows() %>% 
  group_by(SNP) %>%
  arrange(eugene.pvalue) %>% 
  summarize(eugene = paste(glue::glue("{eugene.gene} in {eugene.tissue} (P={eugene.pvalue})"), collapse = ", ") ) %>%
  ungroup()


annotation_gtex <- annotation_raw_gtex %>% 
  bind_rows() %>% 
  group_by(SNP) %>%
  arrange(gtex.pvalue) %>% 
  summarize(gtex = paste(glue::glue("{gtex.gene} in {gtex.tissue} (P={gtex.pvalue})"), collapse = ", ") ) %>%
  ungroup()


## adding annotation to res
## ensembl_gene_stable_id matches with external_gene_name
## renaming gene name
res <-  res %>% 
  left_join(annotation_gene %>% dplyr::select(-ensembl_gene_stable_id, -chr_name), by = c("SNP" = "refsnp_id")) %>%
 # left_join(annotation_eqtlgen, by = c("SNP_tmp" = "SNP")) %>%
  left_join(annotation_gtex, by = c("SNP" = "SNP")) %>%
#  left_join(annotation_eugene , by = c("SNP_tmp" = "SNP")) %>%
  rename(gene = external_gene_name)

## change some gene names

res %<>% mutate(gene = case_when(
  gene == "RP11-141M1.3" ~ "STARD13",
  gene == "AC096570.2" ~ "LINC01830",
  TRUE ~ gene
))


## write out
write_delim(res, path = glue::glue("{DIR_DESCRIPT}/results_table.txt"), delim = " ")



## -------------------------------------------
##   TABLES (XLS)                           --
## -------------------------------------------

#- Descriptive
#- Results (summary)
#- Results (detailed)


## data_summary ----------------------
data_sm_aa <- readr::read_delim(glue::glue("{DIR_DESCRIPT}/descript_table_{NAM_AA}.txt"), delim = " ")
#data_sm_gene <- readr::read_delim(glue::glue("{DIR_DESCRIPT}/descript_table_{NAM_GENE}.txt"), delim = " ")
data_sm_gene1 <- readr::read_delim(glue::glue("{DIR_DESCRIPT}/descript_table_{NAM_GENE_BINARY_1}.txt"), delim = " ")
#data_sm_gene5 <- readr::read_delim(glue::glue("{DIR_DESCRIPT}/descript_table_{NAM_GENE_BINARY_5}.txt"), delim = " ")
data_sm <- rbind(data_sm_gene1[2,], data_sm_aa[2:nrow(data_sm_aa),]) %>% as.data.frame()

## results ---------------------------
## remove rare EBV variants ----------
res <- readr::read_delim(glue::glue("{DIR_DESCRIPT}/results_table.txt"), delim = " ") 


## remove columns with NA only
res %<>% discard(~all(is.na(.x))) %>% dplyr::select(loci, global.outcome, outcome, top_SNP, SNP, Chr, bp, everything())

## remove useless columns
res %<>% dplyr::rename(EAF = Freq) %>%
  select(-SNP_chrpos) %>% 
  as.data.frame() 

## arrange
res %<>% arrange(outcome, p, desc(top_SNP))


## subset
res_sm <- res %>% 
  dplyr::filter(top_SNP | causal_finemap) %>%
  dplyr::select(loci, global.outcome, outcome, SNP, Chr, OR, p, top_SNP, causal_finemap, causal_finemap_prob, effect.allele, EAF, n, gene, consequence_type_tv, gtex)

##dictionary (names(data_sm))
dict_tab1 <- tribble(
  ~`Column Name`, ~`Description`,
  "dataset",  "What type of data",
  "variable",  "What feature is described",
  "n" , "Sample size in that data set", 
  "complete", "complete samples",
  "top_counts", "Frequency of groups (only for categorical vairables)",
  "mean", "Mean of that variable",
  "median", "Median of that variable",
  "sd", "Standard deviation of that variable",
  "min", "Minimum value of that variable",
  "max", "Maximum value of that variable"       
) %>% as.data.frame()
stopifnot(ncol(data_sm) == nrow(dict_tab1))

#names(res)
dict_tab2 <- tribble(
  ~`Column Name`, ~`Description`,
  "loci",  "Locus number",
  "global.outcome",  "EBV dataset",
  "outcome",  "EBV gene/amino acid",
  "top_SNP",  "logical, SNP with locus-wide lowest P-value",
  "causal_finemap",  "logical, indicated as causal by FINEMAP",
  "causal_finemap_prob",  "posterior probability indicated by FINEMAP",
  "SNP",  "SNP identifier",
  "Chr",  "Chromosome",
  "bp",  "Position on build hg19 (37)",
  "effect.allele",  "Coded effect allele",
  "other.allele",  "Coded other allele",
  "EAF",  "Effect allele frequency of SNP",
  "b",  "Effect size in SHCS",
  "se",  "Standard error in SHCS",
  "p",  "P-value in SHCS",
  "outcome.freq",  "Outcome frequency (mean for variant frequency per gene) in SHCS",
  "n",  "Sample size in SHCS",
  "lambda_gc", "Genomic inflation factor",
  "OR",  "Odds ratio (exp(b) for logistic mixed effects model) in SHCS",
  "consequence_type_tv",  "Variant consequence",
  "gene",  "Gene containing SNP",
 # "cis_eqtlgen",  "Gene andissue of cis-eqtl association (P-value of cis-eqtl association), from Vosa et al. 2018",
  "gtex",  "Gene and tissue of eqtl association (P-value of association), from GTEX ",
#  "eugene",  "Gene and tissue of eqtl association (P-value of association), from EUGENE"
) %>% as.data.frame()
stopifnot(ncol(res) == nrow(dict_tab2))

## prep finemap ---------------------------
res_finemap <- fs::dir_ls(DIR_OUTPUT) %>%
  str_subset("results_finemap_") %>%
  str_ignore("validation") %>%
  str_subset(".txt") %>%
  purrr::map_dfr(function(x)
    read_delim(x, delim = " ", col_types = c("dccdddddddcddddccc"))) %>% 
  mutate(global.outcome = case_when(
    global.outcome  %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat" ~ "gene_matrix.max1samp",
    global.outcome  %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_gene_matrix.max5samp.binary.non_synonymous.dat" ~ "gene_matrix.max5samp",
    global.outcome  %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat" ~ "gene_matrix.max0samp",
    global.outcome  %>% str_replace("_debugging", "") %>% str_replace("_validation", "") == "NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat" ~ "aa_variant_matrix"
  )) %>% 
  select(-bp_min, -bp_max, -bp_median, -rank) %>% 
  as.data.frame()

#names(res_finemap)
dict_tab3 <- tribble(
  ~`Column Name`, ~`Description`,
  "Chr", "Chromosome",
  "global.outcome",  "EBV dataset",
  "outcome",  "EBV gene/amino acid",
  "n", "Sample size",
  "bp_from", "Window used for fine-mapping",
  "bp_to", "Window used for fine-mapping",
  "config", "Top ranked causal configuration",
  "prob", "Posterior probability by FINEMAP that config is the causal configuration",
  "log10bf", "log10 Bayes factors",
  "odds",  "Odds of causal configuration",
  "h2", "heritability contribution of SNPs",
  "h2_0.95CI", "95% credible interval of the heritability contribution of SNPs",
  "mean", "joint posterior effect size means",
  "sd", "joint posterior effect size standard deviations"
)%>% as.data.frame()
stopifnot(ncol(res_finemap) == nrow(dict_tab3))

## contigency tables ----------------------

#source("munge/investigate-results-or.R" )
## will return res_contingencytables
## does not need to run, cause not included

## store out to xlsx ----------------------

library(xlsx)
write.xlsx(data_sm, file = glue::glue("{DIR_DESCRIPT}/ms_tables.xlsx"),
           sheetName = "Table 1 (data summary)", append = FALSE, showNA=FALSE, row.names = FALSE) 
write.xlsx(res_sm, file = glue::glue("{DIR_DESCRIPT}/ms_tables.xlsx"),
           sheetName = "Table 2 (results summary)", append = TRUE, showNA=FALSE, row.names = FALSE)

write.xlsx(res, file = glue::glue("{DIR_DESCRIPT}/ms_tables.xlsx"),
           sheetName = "Table S1 (results detailed)", append = TRUE, showNA=FALSE, row.names = FALSE)

write.xlsx(res_finemap, file = glue::glue("{DIR_DESCRIPT}/ms_tables.xlsx"),
           sheetName = "Table S2 (FINEMAP results)", append = TRUE, showNA=FALSE, row.names = FALSE)

#write.xlsx(res_contingencytables, file = glue::glue("{DIR_DESCRIPT}/ms_tables.xlsx"),
#           sheetName = "Table S3 (contigency tables)", append = TRUE, showNA=FALSE, row.names = FALSE)

write.xlsx(dict_tab1, file = glue::glue("{DIR_DESCRIPT}/ms_tables.xlsx"),
           sheetName = "Data dictionary Table 1", append = TRUE, showNA=FALSE, row.names = FALSE)
write.xlsx(dict_tab2, file = glue::glue("{DIR_DESCRIPT}/ms_tables.xlsx"),
           sheetName = "Data dictionary Table 2 and S1", append = TRUE, showNA=FALSE, row.names = FALSE)
write.xlsx(dict_tab3, file = glue::glue("{DIR_DESCRIPT}/ms_tables.xlsx"),
           sheetName = "Data dictionary Table S2", append = TRUE, showNA=FALSE, row.names = FALSE)



## FUMA
#$ scp rueger@grfe3:~/G2G-EBV/data/2_results/lmm_gcta_NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat_BRLF1:p.Glu377Ala.mlma  .
#$ scp rueger@grfe3:~/G2G-EBV/data/2_results/lmm_gcta_NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat_BCRF1.mlma  .
# > http://fuma.ctglab.nl/snp2gene/41127


## -------------------------------------------
##  Results summary
## -------------------------------------------

res <- readr::read_delim(glue::glue("{DIR_DESCRIPT}/results_table.txt"), delim = " ")


res <- fs::dir_ls(DIR_OUTPUT) %>%
  stringr::str_subset("results_lenient_") %>%
  stringr::str_subset(".txt") %>%
  purrr::map_dfr(function(x)
    data.table::fread(x)) %>%
  mutate(p.chr = format(p, digits = 3, scientific = TRUE )) 

#, p.partner.signal.chr = format(p.partner.signal, digits = 3, scientific = TRUE )) %>% filter(p < THRESH)
res <- res %>% mutate(phenotype = case_when(
  str_detect(global.outcome, "gene_matrix.max0samp.depth") ~ "EBV gene variation (fraction)",
  str_detect(global.outcome, "_aa_") ~ "EBV amino acids (binary)",
  str_detect(global.outcome, "gene_matrix.max1samp.binary") ~ "Rare EBV gene variation (binary)",
  str_detect(global.outcome, "gene_matrix.max5samp.binary") ~ "Common EBV gene variation (binary)"
)) %>% dplyr::select(-global.outcome) %>% mutate(OR = exp(b))  %>% filter(phenotype != "gene (frequency)" )


res %>% dplyr::select(outcome, phenotype) %>% 
  distinct() %>% 
  group_by(phenotype) %>% 
  summarize(`#` = n()) %>% 
  kable()

#- Amino acids: 21 phenotypes (24 signals)
#- Amino acids (ns): 8 phenotypes (11 signals)
#- Genes: 1 phenotype
#- Genes (ns): 1 phenotype

qp_summary <-
  ggplot(data = res, aes(outcome.freq, log10(p), colour = phenotype)) + geom_point(alpha = I(0.8)) + geom_point(shape = 1) + xlim(0, NA) + labs(
    title = glue::glue("{nrow(res)} top amino acid associations"),
    subtitle = glue::glue("P < {p.THRESH.lenient}")
  )


## store
png(glue::glue("{DIR_DESCRIPT}/results_summary_maf.png"),
    height = 500, width = 600)
print(qp_summary)
dev.off()




## -------------------------------------------
##   Extensionplot                        --
## -------------------------------------------

## - read all mlma-subset results
## - merge + filter
## - plot all together facet_wrap(outcome), color by type 

## 1. get the top results
res_top <- fs::dir_ls(DIR_OUTPUT) %>%
  str_subset("results_top_") %>%
  str_subset(".txt") %>%
  purrr::map_dfr(function(x)
    fread(x))
res_top_distinct <-
  res_top %>% dplyr::select(outcome, global.outcome) %>% distinct()
OUTCOME <- res_top_distinct$outcome
OUTCOME.GLOBAL <- res_top_distinct$global.outcome


## 2. read all the subset.txt in
window_length <- 1e5

dat <- purrr::map2(OUTCOME, OUTCOME.GLOBAL, function(.x, .y) {
  data_snps_signif <-
    res_top %>% filter(outcome == .x &
                         global.outcome == .y) %>% dplyr::select(Chr, bp)
  
  data <- fs::dir_ls(DIR_SCRATCH) %>%
    str_subset(.x) %>%
    str_ignore("subset.txt") %>%
    str_subset(.y) %>%
    str_ignore("clumped") %>%
    str_ignore("log") %>%
    str_ignore("nosex") %>%
    purrr::map_dfr(function(file)
      fread(file, header = TRUE)  %>% mutate(global.outcome = .y,
                                             outcome = .x)) %>%
    filter(p < 1e-3)
  
  ## filter out clump around significant ones
  data_sig <-
    purrr::map2_dfr(data_snps_signif$Chr, data_snps_signif$bp, function(chr_top, bp_top) {
      data %>%
        filter(Chr == chr_top &
                 bp < (bp_top + window_length) &
                 bp > (bp_top - window_length))
    }) %>% mutate(signif_clump = "signif") %>% distinct()
  
  ## merge back with new column signif_clump, name NAs as "not-signif"
  data <-
    full_join(data, data_sig) %>% mutate(signif_clump = replace_na(signif_clump, "not-signif"))
  
  
  return(data)
  
}) %>%
  bind_rows()


## 3. rename global.outcome
dat <- dat %>% mutate(global.outcome = case_when(
  str_detect(global.outcome, "gene_matrix.max0samp.depth") ~ "EBV gene variation (fraction)",
  str_detect(global.outcome, "_aa_") ~ "EBV amino acids (binary)",
  str_detect(global.outcome, "gene_matrix.max1samp.binary") ~ "Rare EBV gene variation (binary)",
  str_detect(global.outcome, "gene_matrix.max5samp.binary") ~ "Common EBV gene variation (binary)"
))

## 4. Some data transformations for x axis
#dat <- dat %>% 
#  dplyr::arrange(Chr, bp) %>%
#  dplyr::mutate(tmp = 1, bp_cumsum = cumsum(tmp))

library(ggman) ## devtools::install_github("https://github.com/mkanai/ggman")
conv <- ggman:::.convert2posX(dat$Chr, dat$bp, "hg19")

dat$x <- conv$posX

## real distance
#chr_label <-
#  dat %>% group_by(Chr) %>% summarise(median.x = median(bp_cumsum), max.x = max(bp_cumsum))

chr_label <- tibble(labels = conv$labels, breaks = conv$breaks)
  
n_tests <- 458 ## run script in src/eff-number-phenotypes.R for this
THRESH <- p.THRESH.lenient / (n_tests) # 0.05/( 6365738 * n_tests)

## 6. mhtplot with only SIGNIF + labels
peak_labels <- dat %>% filter(signif_clump == "signif") %>% dplyr::select(global.outcome, outcome, x, p, Chr) %>% group_by(global.outcome, outcome, Chr) %>% summarize(x_label = median(x), max.p = max(-log10(p)))

qp_signif <-
  ggplot(data = dat %>% filter(signif_clump == "signif")) +
  geom_hline(yintercept = -log10(THRESH), linetype = 2) +
  geom_point(aes(
    x = x,
    y = -log10(p),
    color = outcome#,
  #  shape = global.outcome
  ),
  alpha = I(0.8)) +
#  scale_color_manual(values = c("black", "black", epfl_red, "black", epfl_turquise, "black")) +  for poster
  scale_color_manual(values = c("black", "black", "black")) + 
  ggrepel::geom_text_repel(
    data = peak_labels,
    aes(x_label, max.p, label = outcome,  color = outcome),
    angle = 45,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    breaks = chr_label$breaks,
    labels = chr_label$labels,
    #,
    expand = c(0.03, 0.03)
  ) +
  xlab("Chromosomes") + 
  guides(color=guide_legend(title="Outcome"), shape=guide_legend(title="Dataset"))

png(glue::glue("{DIR_DESCRIPT}/results_summary_mht_signif.png"), height = 300, width = 900)
print(qp_signif)
dev.off()

ggsave(glue::glue("{DIR_DESCRIPT}/results_summary_mht_signif.pdf"), qp_signif, device = cairo_pdf, height = 5, width = 15)


qp_signif2 <-
  ggplot(data = dat %>% filter(signif_clump == "signif")) +
  geom_hline(yintercept = -log10(THRESH), linetype = 2) +
  geom_point(aes(
    x = x,
    y = -log10(p),
    color = outcome#,
   # shape = global.outcome
  ),
  size = 2,
  alpha = I(0.8)) +
  scale_color_manual(values = c("black", "black", "black")) + 
#  scale_color_manual(values = c("black", "black", epfl_red, "black", epfl_turquise, "black")) + 
  scale_x_continuous(
    breaks = chr_label$breaks,
    labels = chr_label$labels,
    #,
    expand = c(0.03, 0.03)
  ) +
  xlab("Chromosomes") + 
  guides(color=guide_legend(title="Outcome"), shape=guide_legend(title="Dataset"))

ggsave(glue::glue("{DIR_DESCRIPT}/results_summary_mht_signif_nolabel.pdf"), qp_signif2, device = cairo_pdf, height = 5, width = 12)

#c("black", "black", epfl_red, "black", epfl_turquise, "black")
#black             black             red               black             green black
#BALF5             BCRF1             BDLF1             BLLF3             BORF1 BRLF1:p.Glu377Ala 
#4303              7201              7920              5513              6841              8843 

## 5. manhattanplot with shape = signif_clump 


qp_all <- ggplot(data = dat) +
  geom_hline(yintercept = -log10(THRESH), linetype = 2) +
  geom_point(aes(
    x = x,
    y = -log10(p),
    color = outcome,
    shape = signif_clump
  ),
  alpha = I(0.2)) +
  geom_point(
    data = dat %>% filter(signif_clump == "signif"),
    aes(
      x = x,
      y = -log10(p),
      color = outcome,
      shape = signif_clump
    )
  ) +
  facet_wrap( ~ global.outcome) +
  scale_shape_manual(values = c(1, 16)) +
  guides(shape = FALSE) +
  scale_x_continuous(
    breaks = chr_label$breaks,
    labels = chr_label$labels,
    #,
    expand = c(0.03, 0.03)
  ) +
  xlab("Chromosomes")+ 
  guides(color=guide_legend(title="Outcome"))


png(glue::glue("{DIR_DESCRIPT}/results_summary_mht.png"), height = 500, width = 600)
print(qp_all)
dev.off()

ggsave(glue::glue("{DIR_DESCRIPT}/results_summary_mht.pdf"), qp_all, device = cairo_pdf, height = 7, width = 9)



## -------------------------------------------
## print manhattanplot for BDLF1 and BORF1
## for poster
## -------------------------------------------


## only for poster!!!
if (FALSE)
{
  color_highlight <- c(epfl_red, epfl_turquise)
  NAM <-
    "NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat"
  variant <- c("BDLF1", "BORF1")
  
  for (i in 1:length(variant))
  {
    VARIANT <- variant[i]
    out <-
      data.table::fread(glue::glue("{DIR_SCRATCH}/lmm_gcta_{NAM}_{VARIANT}.mlma")) %>% rename(chr = Chr, pos = bp)
    
    ## Lambda QC
    
    lambda_gc <-
      GWAS.utils::genomic_inflation(P = na.omit(out$p))# sel$lambda_gc[i]
    
    
    ##  Q-Q plot -----------
    
    qp.qq <-
      qqplot_wrapper(
        x = out %>% pull(p),
        subtitle = glue::glue("GC lambda = {round(lambda_gc, 2)}"),
        title = glue::glue("{VARIANT} in SHCS")
      )
    
    png(glue::glue(
      "{DIR_OUTPUT}/results_top_qqplot_{NAM}_{VARIANT}_poster.png"
    ))
    print(qp.qq)
    dev.off()
    
    ## MHT plot -----------
    
    png(
      glue::glue(
        "{DIR_OUTPUT}/results_top_manhattan_{NAM}_{VARIANT}_poster.png"
      ),
      width = 900
    )
    qp <-
      mhtplot_wrapper(
        data = out  %>% select(SNP, chr, pos, p),
        title = glue::glue("{VARIANT} GWAS"),
        THRESH1 = THRESH,
        THRESH2 = p.THRESH.lenient,
        color_highlight = color_highlight[i]
      )
    print(qp)
    dev.off()
    
    
    
  }
}

## -------------------------------------------
##   Sample size                   --
## -------------------------------------------



## default sample size
## CHECK in 2_results/*log
#N_GCTA <- 268 
NAM <- c("NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat", 
         "NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat", 
         "NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat", 
         "NC_009334_with_t1_alts_gene_matrix.max5samp.binary.non_synonymous.dat")

N_GCTA <- purrr::map_dfr(NAM, ~ sample_size_extract(DIR = DIR_SCRATCH, phenotype = .x))  

N_GCTA %>% summarize(min(n), max(n), median(n))
# A tibble: 1 x 3
#`min(n)` `max(n)` `median(n)`
#<dbl>    <dbl>       <dbl>
## 1      120      268         264

