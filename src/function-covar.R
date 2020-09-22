## Functions that help with covariates preparations


#' Load covariates file
#'
#' @param file path to file
#'
#' @return data frame
#'
#' @examples 
#' load_prepare_covar(file = glue::glue("{DIR_COVAR}/EBVG2G.covar"))
load_prepare_covar <- function(file = NULL) {
  data <- read_delim(file = file,
                     "\t",
                     escape_double = FALSE,
                     trim_ws = TRUE) %>%
    janitor::clean_names(case = "snake") %>%
    rename_all(tolower)  %>%
    rename(id = iid)
  
  return(data)
}
