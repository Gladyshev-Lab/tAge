#' Get package data file paths
#'
#' This function provides access to example data files included with the tAge package.
#' The data files are located in the inst/extdata directory and include:
#' - Exprs_example.csv: Example expression data
#' - Metadata_example.csv: Example metadata
#' - Gene_list_all_4.6.txt: Gene list for analysis
#'
#' @param filename Character string specifying the name of the data file to retrieve.
#' @return Character string with the full path to the requested data file.
#' @export
#' @examples
#' # Get path to example expression data
#' expr_path <- get_package_data("Exprs_example.csv")
#' 
#' # Get path to example metadata
#' meta_path <- get_package_data("Metadata_example.csv")
#' 
#' # Get path to gene list
#' gene_path <- get_package_data("Gene_list_all_4.6.txt")
get_package_data <- function(filename) {
  if (!file.exists(system.file("extdata", filename, package = "tAge"))) {
    available_files <- list.files(system.file("extdata", package = "tAge"))
    stop("Data file '", filename, "' not found in package. Available files: ", 
         paste(available_files, collapse = ", "))
  }
  return(system.file("extdata", filename, package = "tAge"))
}

#' Load example expression data
#'
#' Loads the example expression data included with the tAge package.
#'
#' @return A data frame with expression data where rows are genes and columns are samples.
#' @export
#' @examples
#' expr_data <- load_example_expression_data()
#' head(expr_data)
load_example_expression_data <- function() {
  expr_path <- get_package_data("Exprs_example.csv")
  return(read.csv(expr_path, row.names = 1))
}

#' Load example metadata
#'
#' Loads the example metadata included with the tAge package.
#'
#' @return A data frame with sample metadata.
#' @export
#' @examples
#' meta_data <- load_example_metadata()
#' head(meta_data)
load_example_metadata <- function() {
  meta_path <- get_package_data("Metadata_example.csv")
  return(read.csv(meta_path, row.names = 1))
}

#' Load gene list
#'
#' Loads the gene list included with the tAge package.
#'
#' @return A character vector of gene identifiers.
#' @export
#' @examples
#' gene_list <- load_gene_list()
#' head(gene_list)
load_gene_list <- function() {
  gene_path <- get_package_data("Gene_list_all_4.6.txt")
  gene_list <- readLines(gene_path)
  gene_list <- gene_list[gene_list != ""]
  gene_list <- trimws(gene_list)
  return(gene_list)
}
