# Build function to check info
#
# Build function based on some argument
#
# @param keytype keytype in annotation package
#
# @return return a function object to check symbol whether in specific type
#
#
is_annotation <- function(keytype) {
  .is_annotation <- function(genes, species = c('Hs', 'Mm')) {
    species <- match.arg(species)
    if (species == "Hs") {
      org <- Homo.sapiens::Homo.sapiens
    } else {
      org <- Mus.musculus::Mus.musculus
    }
    universe <- AnnotationDbi::keys(org, keytype = keytype)
    genes %in% universe
  }
}

#' Check Gene id
#'
#' Check Gene id whether is matched with specific id type
#'
#' @param genes gene id to check
#' @param species The species of gene. only support Hs or Mm
#'
#' @aliases is_symbol is_entrez is_ensembl is_refseq is_accnum is_ucsc
#'
#' @return return a logical vector
#'
#' @name check_id
#' @rdname check_id
#'
#' @examples
#' \dontrun{
#'
#' is_symbol('TP53', species = 'Hs')
#' is_symbol(c('TP53', 'SJDSJ'))
#'
#' }
#'
NULL


#' @rdname check_id
#'
#' @export
#'
is_symbol <- is_annotation("SYMBOL")

#' @rdname check_id
#'
#' @export
#'
is_entrez <- is_annotation("ENTREZID")

#' @rdname check_id
#'
#' @export
#'
is_ensembl <- is_annotation("ENSEMBL")

#' @rdname check_id
#'
#' @export
#'
is_refseq <- is_annotation("REFSEQ")

#' @rdname check_id
#'
#' @export
#'
is_accnum <- is_annotation("ACCNUM")

#' @rdname check_id
#'
#' @export
#'
is_ucsc <- is_annotation("UCSCKG")


#' Check gene id type
#'
#' Check gene id type by specific gene
#'
#'
#' @param symbol The gene symbol to check
#' @param species The species of gene. only support Hs or Mm
#' @param aliasesBeforeSymbols check aliases before Symbols. If TRUE, when aliases is equal to Symbol, will return 'aliases'. Default FALSE.
#' @param priority The priority to check gene id type. you can set a vector of id type.
#'
#' @seealso
#' \code{\link[AnnotationDbi]{ACCNUM}}
#' \code{\link[Homo.sapiens]{Homo.sapiens}}
#' \code{\link[Mus.musculus]{Mus.musculus}}
#'
#' @return return a vector with symbol type
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' get_geneid_type('TP53', species = 'Hs')
#' get_geneid_type(c('TP53', 'ADSL'), species = 'Hs')
#'
#' }
#'
#'
get_geneid_type <- function(symbol, species = c('Hs', 'Mm'), aliasesBeforeSymbols = FALSE, priority = NULL) {
  keytypes <- c("SYMBOL", "ALIAS", "ENTREZID", "ENSEMBL", "REFSEQ", "ACCNUM", "UCSCKG",
                "TXID", "ENSEMBLTRANS", "ENSEMBLPROT", "UNIPROT", "UNIGENE", "EXONID", "ENZYME")

  if (aliasesBeforeSymbols) {
    if (!is.null(priority) && priority != "SYMBOL") {
      keytypes[1:2] <- keytypes[2:1]
    } else {
      message("You set priority to 'SYMBOL' so aliasesBeforeSymbols will be ignored.")
    }
  }

  if (!is.null(priority)) {
    if (!priority %in% keytypes) {
      stop(priority, " does not exist")
    } else {
      keytypes <- c(priority, setdiff(keytypes, priority))
    }
  }

  .call <- function(keytype) {
    if (any(is.na(nams))) {
      nams[is.na(nams) & is_annotation(keytype)(symbol, species)] <<- keytype
    }
  }

  nams <- rep(NA, length(symbol))
  sapply(keytypes, .call)
  setNames(nams, symbol)
}


