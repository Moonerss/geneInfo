
# function midoficated based on TCGAbiolinks package

#' @title Get gene information from biomaRt
#' @description Get  information from biomaRt
#'
#' @param species choose the species, 'Hs' or 'Mm'
#' @param hg19 It is specific to \code{species = 'Hs'}, whether download hg109 information
#' @param mm10 It is specific to \code{species = 'Mm'}, whether download mm10 information
#' @param ... Other argument of \code{\link[biomaRt]{useEnsembl}}
#'
#' @return return a data.frame object of gene coordinate information
#'
#' @seealso https://www.ensembl.org/info/website/archives/assembly.html
#'
#' @importFrom biomaRt getBM useEnsembl listDatasets
#'
#' @export
#'
get_gene_coord_from_ensembl <- function(species = c('Hs', 'Mm'), hg19 = FALSE, mm10 = FALSE, ...) {
  species <- match.arg(species)

  ## check version info
  dot_list <- list(...)
  if ('version' %in% names(dot_list)) {
    if (species == 'Mm') {
      if (is.null(dot_list$version)) {
        if (isTRUE(mm10)) {
          version <- '102'
        } else {
          version <- dot_list$version
        }
      } else {
        version <- dot_list$version
      }
    } else {
      version <- dot_list$version
    }
  } else {
    if (species == 'Mm') {
      if (isTRUE(mm10)) {
        version <- '102'
      } else {
        version <- NULL
      }
    } else {
      version <- NULL
    }
  }

  tries <- 0L
  msg <- character()
  ## try 3 times to get info
  while (tries < 3L) {
    gene.location <- tryCatch({
      ## ensembl info
      if (species == 'Hs') {
        dataset <- 'hsapiens_gene_ensembl'
        if (isTRUE(hg19)) {
          host <- "https://grch37.ensembl.org"
        } else {
          host <- "www.ensembl.org"
        }
      } else {
        dataset <- 'mmusculus_gene_ensembl'
        if (isTRUE(mm10)) {
          host <- 'https://nov2020.archive.ensembl.org'
        } else {
          host <- "www.ensembl.org"
        }
      }

      mirror <- list(NULL, "useast", "uswest", "asia")[[tries + 1]]
      ensembl <- tryCatch({
        message(ifelse(is.null(mirror),
                       paste0("Accessing ", host, " to get gene information"),
                       paste0("Accessing ", host," (mirror ", mirror,")")))
        useEnsembl("ensembl", dataset = dataset, host = host, mirror = mirror, version = version, ...)
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





#' @title Get gene information from bioconductor
#' @description Get  information from bioconductor \code{db} package and corresponding \code{Txdb} package
#'
#' @param species choose the species, 'Hs' or 'Mm'
#' @param hg19 It is specific to \code{species = 'Hs'}, whether download hg109 information
#' @param mm10 It is specific to \code{species = 'Mm'}, whether download mm10 information
#'
#' @return return a data.frame object of gene coordinate information
#'
#' @seealso
#' \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}
#' \code{\link[TxDb.Hsapiens.UCSC.hg19.knownGene]{TxDb.Hsapiens.UCSC.hg19.knownGene}}
#' \code{\link[TxDb.Hsapiens.UCSC.hg38.knownGene]{TxDb.Hsapiens.UCSC.hg38.knownGene}}
#' \code{\link[org.Mm.eg.db]{org.Mm.eg.db}}
#' \code{\link[TxDb.Mmusculus.UCSC.mm10.knownGene]{TxDb.Mmusculus.UCSC.mm10.knownGene}}
#' \code{\link[TxDb.Mmusculus.UCSC.mm39.knownGene]{TxDb.Mmusculus.UCSC.mm39.knownGene}}
#'
#' @importFrom stats aggregate
#' @importFrom dplyr left_join
#' @importFrom rlang .data
#'
#' @export
#'
get_gene_coord_from_txdb <- function(species = c('Hs', 'Mm'), hg19 = FALSE, mm10 = FALSE) {

  species <- match.arg(species)

  if (species == 'Hs') {
    message('get gene annotation from: org.Hs.eg.db')
    orgdb <- org.Hs.eg.db::org.Hs.eg.db
    if (isTRUE(hg19)) {
      message('get gene coordinate from: TxDb.Hsapiens.UCSC.hg19.knownGene')
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    } else {
      message('get gene coordinate from: TxDb.Hsapiens.UCSC.hg38.knownGene')
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }
  } else if (species == 'Mm') {
    message('get gene annotation from: org.Mm.eg.db')
    orgdb <- org.Mm.eg.db::org.Mm.eg.db
    if (isTRUE(mm10)) {
      message('get gene coordinate from: TxDb.Mmusculus.UCSC.mm10.knownGene')
      txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    } else {
      message('get gene coordinate from: TxDb.Mmusculus.UCSC.mm39.knownGene')
      txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene::TxDb.Mmusculus.UCSC.mm39.knownGene
    }
  }
  canonical.chr <- paste0("chr", c(1:22, "X", "Y"))
  symbols <- AnnotationDbi::keys(orgdb, keytype = "ENTREZID")
  symbol_output <- AnnotationDbi::select(orgdb,
                                         keys = symbols,
                                         keytype = "ENTREZID",
                                         columns = c('SYMBOL', 'ENSEMBL', 'GENETYPE', 'ALIAS'))
  location_output <- as.data.frame(GenomicFeatures::genes(txdb))

  gene_location <- symbol_output %>%
    dplyr::left_join(location_output, by = c('ENTREZID' = 'gene_id')) %>%
    dplyr::filter(!is.na(.data$start))

  paste_alias <- aggregate(ALIAS ~ ENTREZID, data = gene_location[, c('ENTREZID', 'ALIAS')], FUN = function(x) paste(x, collapse = ","))

  final_gene_location <- gene_location[!duplicated(gene_location$ENTREZID), ]
  final_gene_location <- final_gene_location %>%
    dplyr::filter(.data$seqnames %in% canonical.chr) %>%
    dplyr::select(-.data$ALIAS) %>%
    left_join(paste_alias, by = 'ENTREZID') %>%
    dplyr::select(.data$SYMBOL, .data$seqnames, .data$start, .data$end, .data$strand, .data$ENSEMBL, .data$ENTREZID, .data$ALIAS, .data$GENETYPE)

  colnames(final_gene_location) <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position",
                                     "strand", "ensembl_gene_id", "entrezgene_id", "external_gene_name",
                                     'gene_biotype')

  return(final_gene_location)
}
