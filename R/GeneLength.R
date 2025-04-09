#' @title Get Gene exon length information
#' @description
#' Get gene info from \code{OrganismDb} package
#'
#' @param symbols character vector of HGNC gene symbols
#' @param species The species of gene. only support Hs or Mm
#'
#' @note
#' Only consider the gene on chr 1-22, chrX and chrY
#'
#' @export
#' @rdname get_exons_info
#' @examples
#' \dontrun{
#' get_exons_info('TP53')
#' get_exons_info('Pzp', species = 'Mm')
#' }
#'
get_exons_info <- function(symbols = NULL, species = c("Hs", "Mm")) {
  species <- match.arg(species)
  if (species == "Hs") {
    orgdb <- Homo.sapiens::Homo.sapiens
  } else {
    orgdb <- Mus.musculus::Mus.musculus
  }

  canonical.chr <- paste0("chr", c(1:22, "X", "Y"))
  if (is.null(symbols)) {
    symbols <- AnnotationDbi::keys(orgdb, keytype = "SYMBOL")
  }
  out <- AnnotationDbi::select(orgdb,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("EXONID", "EXONRANK", "EXONCHROM", "EXONSTART", "EXONEND")
  ) %>%
    dplyr::filter(EXONCHROM %in% canonical.chr) %>%
    dplyr::group_by(EXONCHROM, SYMBOL) %>%
    dplyr::summarize(
      n_allowed_exons = length(unique(EXONRANK)),
      n_total_exons = length(unique(EXONID)),
      exon_start = min(EXONSTART),
      exon_end = max(EXONEND),
      gene_length = exon_end + 1 - exon_start
    ) %>%
    dplyr::rename(gene = SYMBOL, exon_chr = EXONCHROM)

  return(out)
}

#' @title Get allowedExons number
#' @description
#' Get allowedExons from \code{OrganismDb} package
#'
#' @inheritParams get_exons_info
#'
#' @note
#' Only consider the gene on chr 1-22, chrX and chrY
#'
#' @importFrom stats setNames
#' @export
#' @examples
#' \dontrun{
#' allowedExons('TP53')
#' allowedExons('Pzp', species = 'Mm')
#' }
#'
allowedExons <- function(symbols = NULL, species = c("Hs", "Mm")) {
  out <- get_exons_info(symbols = symbols, species = species)
  stats::setNames(out$n_allowed_exons, out$gene)
}


#' @title Get totalExons number
#' @description
#' Get totalExons from \code{OrganismDb} package
#'
#' @inheritParams get_exons_info
#'
#' @note
#' Only consider the gene on chr 1-22, chrX and chrY
#'
#' @importFrom stats setNames
#' @export
#' @examples
#' \dontrun{
#' totalExons('TP53')
#' totalExons('Pzp', species = 'Mm')
#' }
#'
totalExons <- function(symbols = NULL, species = c("Hs", "Mm")) {
  out <- get_exons_info(symbols = symbols, species = species)
  stats::setNames(out$n_total_exons, out$gene)
}

#' @title Get gene length
#' @description
#' Get gene length from \code{OrganismDb} package
#'
#' @inheritParams get_exons_info
#'
#' @note
#' Only consider the gene on chr 1-22, chrX and chrY
#'
#' @importFrom stats setNames
#' @export
#' @examples
#' \dontrun{
#' geneLength('TP53')
#' geneLength('Pzp', species = 'Mm')
#' }
#'
geneLength <- function(symbols = NULL, species = c("Hs", "Mm")) {
  out <- get_exons_info(symbols = symbols, species = species)
  stats::setNames(out$gene_length, out$gene)
}
