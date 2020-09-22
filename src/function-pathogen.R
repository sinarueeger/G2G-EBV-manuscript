




## -- Create fraction of counts when gene --------------------------------------------------------

#' Title
#'
#' @param data dataset of phentype, should have one column names id (wiht the identifiers) and the remaining columns gene counts
#' @param path_gene_length path to gene length
#'
#' @return dataset of same size as data
#'
#' @examples
#' data <- fread(glue::glue("{DIR_PATHOGEN}/NC_009334_with_t1_alts_gene_matrix.max0samp.depth_corr_counts.non_synonymous.dat"))%>% rename(id = sample)
#' gene_counts2fraction(data = data, path_gene_length = glue::glue("{DIR_PATHOGEN}/NC_009334_with_t1_alts_gene_summary.max0samp.dat"))

gene_counts2fraction <-
  function(data = NULL,
           path_gene_length = NULL)
  {
    gene_length <- fread(path_gene_length)
    gene_long <- gather(data, gene, count,-id)
    
    gene_long_rel <-
      gene_long %>% left_join(gene_length %>% select(gene, total, synonymous), by = "gene") %>%
      dplyr::mutate(nonsynonymous = total - synonymous) %>%
      dplyr::mutate(rel = count / nonsynonymous)
    
    gene_rel <-
      spread(gene_long_rel %>% select(id, gene, rel), gene, rel)
    
    if (!identical(dim(gene_rel), dim(data)))
      stop("dimension of output data is not the same as input data.")
    
    return(gene_rel)
  }


## -- Load pathogen data --------------------------------------------------------

#' Title
#'
#' @param files_coverage
#' @param files_data
#' @param debugging
#' @param path_gene_length path to gene length file (defined in setup.R)
#' @param counts2frac logic, for genes, should binary variables be transformed into fraction? 
#' @param return_count logic, return count data in the case of gene stuff
#'
#' @return
#' @export
#'
#' @examples
load_prepare_pathogen <-
  function(files_coverage = NULL,
           files_data = NULL,
           debugging = FALSE,
           path_gene_length = NULL,
           counts2frac = FALSE,
           return_count = FALSE) {
    ##  files_coverage, path to file with coverage information
    ##  files_data, path to file with aminoacid, gene, gen information
    
    ## the structure fo the data is always
    ## sample  somestuff
    
    #  browser()
    
    if (!is.na(files_coverage))
    {
      ## coverage
      dat_coverage <- fread(files_coverage, na.strings = "-9.0") %>%
        janitor::clean_names(case = "snake") %>%
        rename_all(tolower) %>%
        rename(depth_above_10x = depth_above_10) %>%
        rename(id = sample) %>%
        select(id, depth_above_10x)
      
      ## depth_above_6x = depth_above6x,  we dont care about this
    }
    
    ## data
    dat_aa <- data.table::fread(files_data, na.strings = "-9.0") %>%
      rename(id = sample) #%>% janitor::clean_names()
    ## If there is less than 6 reads mapping onsite, value is -9. If there has been an unexpected error, value is -1.
    ## I recommend to convert -9 to NAs.
    
    ## check if any -1
    if (any(dat_aa < 0, na.rm = TRUE))
      stop("negative values in input file")
    
    ## debugging
    if (debugging) {
      set.seed(3)
      dat_aa <- dat_aa[, c(1, sample(2:ncol(dat_aa), 15))]
    }
    
    outcomes <- names(dat_aa)[-1]
    
    ## if genecounts, turn into relative gene counts
    ## if gene, calculate relative proportion of counts
    if (counts2frac)
    {
      dat_counts <- dat_aa
      dat_aa <-
        gene_counts2fraction(data = dat_counts, path_gene_length = path_gene_length)
    }
    
    
    ## no need to add info about batch
    if (!is.na(files_coverage))
    {
      ## add info about batch to dat_aa
      dat_not_joined <-
        anti_join(dat_coverage, dat_aa, by = c("id"))
      dat_aa_coverage <-
        dat_aa %>% left_join(dat_coverage, by = c("id"))
    } else{
      dat_aa_coverage <- dat_aa
    }
    
    if (return_count)
    {
      return(list(
        data = dat_aa_coverage,
        outcome = outcomes,
        data_counts = dat_counts
      ))
    } else{
      return(list(data = dat_aa_coverage, outcome = outcomes))
    }
  }



## -- Quality Control ------------------------------------------------------------------------

#' Title
#'
#' @param data_wide
#' @param cut
#'
#' @return
#' @export
#'
#' @examples
apply_qc_pathogen <- function(data_wide, cut = 0.19) {
  ## remove individuals
  ## cut: apparently not really cut at 0.2, but rather at 0.19
  ## only filter for EUR individuals, no QC info available for Sanger individuals
  
  if ("depth_above_10x" %in% names(data_wide$data))
  {
    data_wide_qc <-
      data_wide$data %>% dplyr::filter(depth_above_10x >= cut)
  } else{
    data_wide_qc <- data_wide$data ## we dont care about "above_10"
  }
  
  ## remove individuals with only NAs
  NA_rowwise <-
    apply(data_wide_qc[, -1], 1, function(x)
      sum(is.na(x))) < (ncol(data_wide_qc) - 1)
  data_wide_qc <- data_wide_qc[NA_rowwise, ]
  
  ## remove rare variants?? how to do generally?
  
  return(list(data = data_wide_qc, outcome = data_wide$outcome))
  
}









## -- Calculate PCs ------------------------------------------------------------------------


#' Extract PCs
#' This is a function to compute the phylogenetic PCs for the HBV project
#' Guide on constructing pPCs: https://gist.github.com/sinarueeger/718c7ec3cee473ad0beed4feea4d0e24
#'
#' @param path_tree path to tree
#' @param ppca logical, returning phylogenetic PCs
#' @param pca logical, returning PCs
#' @param data_trait dataset with id and variants as columns
#' @param debugging logical
#' @param validation logical
#' @param n.pc number of PCs returned
#' @param loadings logical, return loadings ("rotation") instead of pcs
#' @param filter_threshold maf filter, everything < filter_threshold will be removed
#'
#' @return dataframe with ID and PCs as columns
#' @examples
#' pathogen_pca(
#' path_tree = glue::glue("{DIR_PATHOGEN}/../EBNA_3A.hiv.msa.manually_trimmed.fasta.treefile") # ape::read.tree(path_tree_EBNA_3A),
#' data_trait = data_pathogen,
#' ppca = TRUE
#' )


pathogen_pca <- function(path_tree = NULL,
                         data_trait,
                         ppca = FALSE,
                         pca = TRUE,
                         n.pc = 4,
                         debugging = FALSE,
                         validation = FALSE,
                         loadings = FALSE, 
                         eigenvalues = FALSE, 
                         filter_threshold = 0.05)
  
{
  if (ppca)
  {
    ## 1) TREE --------------------------------------
    ## ----------------------------------------------
    tree <- ape::read.tree(path_tree)
    #plot(tree)
    
    ## check if matrix is singular
    #solve(ape::vcv.phylo(tree))
    
    ## 2) Y matrix ----------------------------------
    ## ----------------------------------------------
    
    tree_subset <-
      ape::keep.tip(tree, as.character(data_trait$id))
    #    plot(tree_subset)
    
    ## reorder rows, same tiplabels
    id_tree <- tibble(id = as.numeric(tree_subset$tip.label))
    
    ## turn data into matrix and reorder rows according to tips
    Y <- data_trait %>% right_join(id_tree)
    
    stopifnot(identical((tree_subset$tip.label), as.character(Y$id)))
    
  } else {
    Y <- data_trait
  }
  
  id_pathogen <- Y$id
  
  ## remove ID
  Y_no_id <- Y %>% select(-id) %>% as.matrix()
  
  ## remove columns with only wild type
  Y_no_id <- Y_no_id[,which(apply(Y_no_id, 2, function(x) mean(x, na.rm=T)) != 0)]
  
  ## Remove variant between 0.05 and 0.95 allele frequency
  Y_no_id.freq <- parallel::mclapply(1:ncol(Y_no_id), function(x){mean(Y_no_id[,x], na.rm = T)}, mc.cores = 20) %>% unlist()
  Y_no_id <- Y_no_id[,which(Y_no_id.freq > filter_threshold)]
  #>>> only 838
  
  ## only non varying ones
  Y_no_id.var <- parallel::mclapply(1:ncol(Y_no_id), function(x){var(Y_no_id[,x], na.rm = T)}, mc.cores = 20) %>% unlist()
  Y_no_id <- Y_no_id[, which(Y_no_id.var != 0)]
  
  ## impute 
  Y_no_id4logitPCA <- round(missMDA::imputePCA(Y_no_id)$completeObs)
  
  ## initiate out
  out <- data.frame(id = id_pathogen)
  
  ## 4) as comparison: PCA ------------------------
  ## ----------------------------------------------
  if (pca)
  {
  
    if (debugging) {
      n.pc <- 2
    }
    
    ## FactoMineR::PCA ---------------------------
    
    pca_logistic_raw <-
      logisticPCA::convexLogisticPCA(Y_no_id4logitPCA, k = n.pc)
    ## PCs (PC scores)
    ## U = loadings
    
    
    if (loadings) {
      ## return loadings
      pca_out_loadings <-
        pca_logistic_raw$U %>% as.data.frame()  %>% select(1:n.pc)
      out <- pca_out_loadings
    } else {
      if (eigenvalues) {
    
        ## only for eigenvalues, a bit dodgy tho to compute this when logistic PCA was used
        pca_raw <-
          FactoMineR::PCA(Y_no_id4logitPCA, scale.unit = FALSE, graph = FALSE)
        
        ## return eigenvalues
        pca_out_eigenvalues <- (pca_raw$eig)#[1:n.pc]
        out <- pca_out_eigenvalues
        
      } else{
      ## default, return pcs
      pca_out <- pca_logistic_raw$PCs %>% as.data.frame() %>% select(1:n.pc)
      names(pca_out) <- paste0("PC", 1:n.pc)
      out <- cbind(out, pca_out)
      
      }
    }
    
  }
  
  
  
  
  
  ## 3) pPCA --------------------------------------
  ## ----------------------------------------------
  
  
  if (ppca)
  {
    ## from Nimisha
    vir_4d <- phylobase::phylo4d(tree_subset, Y_no_id)
    vir_pca <- adephylo::ppca(
      vir_4d,
      scale = FALSE,
      scannf = FALSE,
      nfposi = 16,
      method = "Abouheif"
    )
    
    
    ppca_out <- vir_pca$li %>% select(1:n.pc)
    
    # plot(vir.pca$li[,1:2])
    #  plot(vir.pca$ls[,1:2])
    
    names(ppca_out) <- paste0("pPC", 1:n.pc)
    
    out <- cbind(out, ppca_out)
  }
  
  
  
  
  return(out)
  
  
}









## -- Visualise data ------------------------------------------------------------------------


#' Title
#'
#' @param data
#' @param filename
#' @param path_out
#' @param format
#'
#' @return
#' @export
#'
#' @examples
summarise_pathogen <-
  function(data,
           filename = NULL,
           path_out = "dontuseme",
           format = "png") {
    ##filename of pathogen
    ## pathout of results
    
    ## if format == ""pdf is used, no stegasaur can be used
    
    ## turn into binary factors if needed
    length.unique <-
      apply(data$data[, data$outcome], 2, function(x)
        length(unique(na.omit(x))))
    type_outcome <- case_when(
      any(length.unique > 2) ~ "continuous",
      all(length.unique <= 2) ~ "binary",
      TRUE ~ "no variation"
    )
    
   # print(type_outcome)
    if (type_outcome == "binary") {
      ## write out summary
      sm <-
        apply(data$data %>% select(data$outcome), 2, function(x)
          sum(x == 1, na.rm = TRUE) / length(x)) %>% t() %>% as_tibble() %>% gather()
      #  write(sm, file = "tnmp.txt")
      
      
      ## histogram of frequencies
      
      qp <- ggplot(data = sm) + geom_histogram(aes(x = value)) +
        xlab("Frequency of 1 per variant") +
        labs(
          title = filename,
          subtitle = "(binary variable)",
          caption = glue::glue("{nrow(sm)} variables")
        )
      
      
    }
    
    if (type_outcome == "continuous") {
      ## write out summary
      sm <-
        apply(data$data %>% select(data$outcome), 2, function(x)
          mean(x, na.rm = TRUE)) %>% t() %>% as_tibble() %>% gather()
      
      
      ## histogram of mean
      qp <- ggplot(data = sm) + geom_histogram(aes(x = value)) +
        xlab("Mean per SNP") +
        labs(
          title = filename,
          subtitle = "(continuous variable)",
          caption = glue::glue("{nrow(sm)} variables")
        )
    }
    
    if (format == "df") {
      return(sm)
    }
    
    if (format == "png") {
      ## replacing dots with underscores
      path_out <-
        paste0(str_replace_all(path_out, "[\\.]", "\\_"), ".", format)
      
      png(path_out)
      print(qp)
      dev.off()
      
      ## now encode the summary
      path_out <-
        paste0(str_replace_all(path_out, "[\\.]", "\\_"), ".txt")
      write_delim(sm, path_out)
    }
  }
