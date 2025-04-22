
# function taken from TCGAbiolinks package

#' @title Get gene information from biomaRt
#' @description Get  information from biomaRt
#'
#' @param species choose the species, 'Hs' or 'Mm'
#' @param hg19 It is specific to \code{species = 'Hs'}, whether download hg109 information
#' @param ... Other argument of \code{\link[biomaRt]{useEnsembl}}
#'
#' @return return a data.frame object of gene coordinate information
#'
#' @importFrom biomaRt getBM useEnsembl listDatasets
#' @export
#'
get_gene_coord_from_ensembl <- function(species = c('Hs', 'Mm'), hg19 = FALSE, ...) {
  species <- match.arg(species)

  tries <- 0L
  msg <- character()
  ## try 3 times to get info
  while (tries < 3L) {
    gene.location <- tryCatch({
      ## ensembl info
      if (species == 'Hs') {
        dataset <- 'hsapiens_gene_ensembl'
        if (isTRUE(hg19)) {
          host <- "grch37.ensembl.org"
        } else {
          host <- "www.ensembl.org"
        }
      } else {
        dataset <- 'mmusculus_gene_ensembl'
        host <- "www.ensembl.org"
      }

      mirror <- list(NULL, "useast", "uswest", "asia")[[tries + 1]]
      ensembl <- tryCatch({
        message(ifelse(is.null(mirror),
                       paste0("Accessing ", host, " to get gene information"),
                       paste0("Accessing ", host," (mirror ", mirror,")")))
        useEnsembl("ensembl", dataset = dataset, host = host, mirror = mirror, ...)
      }, error = function(e) {
        message(e)
        return(NULL)
      })

      attributes <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position",
                      'band', "strand", "ensembl_gene_id", "entrezgene_id", "external_gene_name",
                      'gene_biotype')

      db.datasets <- listDatasets(ensembl)
      if (species == 'Hs') {
        description <- db.datasets[db.datasets$dataset == "hsapiens_gene_ensembl",]$description
      } else {
        description <- db.datasets[db.datasets$dataset == "mmusculus_gene_ensembl",]$description
      }

      message(paste0("Downloading genome information (try:", tries,") Using: ", description))

      chrom <- c(1:22, "X", "Y")
      gene.location <- getBM(attributes = attributes, filters = c("chromosome_name"), values = list(chrom), mart = ensembl)

      gene.location
    }, error = function(e) {
      msg <<- conditionMessage(e)
      tries <<- tries + 1L
      NULL
    })
    if(!is.null(gene.location)) break
  }
  if (tries == 3L) stop("failed to get URL after 3 tries:", "\n  error: ", msg)
  ## modified info
  gene.location$strand[gene.location$strand == 1] <- "+"
  gene.location$strand[gene.location$strand == -1] <- "-"
  gene.location$chromosome_name <- paste0("chr",gene.location$chromosome_name)
  return(gene.location)
}


