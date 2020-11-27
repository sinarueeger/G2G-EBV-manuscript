## calculate genomic lambdas for all mlma outputs
## --------------------------------------------------


##////////////////////////////////////////////////////////////////
##                           0. Setup                           //
##////////////////////////////////////////////////////////////////


args <- commandArgs(trailingOnly = TRUE)


## some default values
PATHOGEN <- args[1] ## 1:length(NAMS)



## This R script creates all the variables
## takes counter_for_pathogen_outcome as input!
DIR_SRC <- here::here("src")
source(glue::glue("{DIR_SRC}/setup.R"))



## read in plan and add sample size
plan <-
  read_delim(glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"), delim = " ") %>% filter(include.gcta)# %>% slice(1:3)



## ////////////////////////////////////////////////////////////
## define two sets of SNPs: 
## - lenient QC (what we have)
## - strict QC (what we do below)
## ////////////////////////////////////////////////////////////
 
## done in make_2_results.sh (src/post_QC)

## Define HOST data
HOST <-
  plan %>% select(host) %>% distinct() %>% assertr::verify(nrow(.) == 1) %>% pull(host)
NAM_HOST <- fs::path_file(HOST)##  glue::glue("{(DIR_HOST)}/"))

## read the fam file (for later purposes)
bim.file <-
  read_delim(glue::glue("{DIR_SCRATCH}/{NAM_HOST}_QC.bim"),
             delim = "\t",
             col_names = FALSE)



##////////////////////////////////////////////////////////////////
##                  2. GENOMIC LAMBDA                           //
##////////////////////////////////////////////////////////////////

## define smallest set (unique outcome and global outcome).
## E.g. if there are multiple loci present within one outcome, 
## only one GC analysis or MHT plot is produced. 

sel <- plan %>% dplyr::select(global.outcome, outcome, outcome.freq) %>%
  unique() %>% mutate(lambda_gc = NA, lambda_gc_sub = NA, sig_lambda_gc = FALSE)


## Loop through all possible results and compute genomic lambda

sel_raw <- parallel::mclapply(1:nrow(sel), function(i){
  
  NAM <- sel[i, "global.outcome"]
  VARIANT <- sel[i, "outcome"]
  FREQ <- sel[i, "outcome.freq"]
  
  gc <-  tryCatch(
      expr = {
        out <- data.table::fread(glue::glue("{DIR_SCRATCH}/lmm_gcta_{NAM}_{VARIANT}.mlma")) %>% rename(chr = Chr, pos = bp)
        
        ## calculate GC lambda
        lambda_gc <- GWAS.utils::genomic_inflation(P = na.omit(out$p))
        sel$lambda_gc[i] <- lambda_gc
        
        ## subset of SNPs
        
        lambda_gc_sub <- GWAS.utils::genomic_inflation(P = na.omit(out$p[out$SNP %in% bim.file$X2]))
        sel$lambda_gc_sub[i] <- lambda_gc_sub
        
        sel[i,]
        
        
      },
      error = function(e){
        message('Caught an error!')
        print(NULL)
      },
      warning = function(w){
        message('Caught an warning!')
        print(NULL)
      },
      finally = {
        message('All done, quitting.')
      }
    )    
  
  return(gc)
  
}, mc.cores = 20)

sel_lambda <- sel_raw %>% dplyr::bind_rows() 

NAM <- sel[1, "global.outcome"]
save(sel_lambda, file = glue::glue("{DIR_SCRATCH}/genomic_lambda_{NAM}.RData"))
range(sel_lambda$lambda_gc, na.rm = TRUE)




## EVALUATION
## -------------

### 
load(glue::glue("{DIR_SCRATCH}/genomic_lambda_NC_009334_with_t1_alts_gene_matrix.max1samp.binary.non_synonymous.dat.RData"))
gc_gene1 <- sel_lambda

load(glue::glue("{DIR_SCRATCH}/genomic_lambda_NC_009334_with_t1_alts_aa_variant_matrix.non_synonymous.binary.dat.RData"))
gc_aa <- sel_lambda

gc_wide <- rbind(gc_gene1, gc_aa)
gc <- gc_wide %>% tidyr::pivot_longer(cols = starts_with("lambda_gc"), names_to = "qc_done", values_to = "lambda_gc")

## range

gc %>% group_by(qc_done) %>% summarize(mean = print(format(mean(lambda_gc), dig = 2)), 
                                       min = print(format(min(lambda_gc), dig = 2)), 
                                       max = print(format(max(lambda_gc), dig = 3)))

## plotting the same
ggplot(data = gc) + geom_violin(aes(y = lambda_gc, x = global.outcome, group = interaction(global.outcome, qc_done), color = qc_done)) + coord_flip()
ggplot(data = gc) + geom_boxplot(aes(y = lambda_gc, x = global.outcome, group = interaction(global.outcome, qc_done), color = qc_done)) + coord_flip()

ggplot(data = gc_wide, aes(x = lambda_gc, y = (lambda_gc - lambda_gc_sub), color = global.outcome)) + geom_point() + 
  geom_abline(intercept = -1, slope = 1)


## outcome freq 
ggplot(data = gc, aes(x = outcome.freq, y = lambda_gc, color = global.outcome)) + geom_point(aes(shape = qc_done)) + scale_shape_manual(values = c(16, 1))
