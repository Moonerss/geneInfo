# convert gene id
# convert gene id between different type
# @param x The first gene type
# @param y The gene type to convert
#
# @return return a function can do gene id convert
#
convert_id <- function(x, y) {
  .convert_id <- function(keys = NULL, species = c("Hs", "Mm"), multiVals = "first") {
    species <- match.arg(species)
    if (species == "Hs") {
      org <- Homo.sapiens::Homo.sapiens
    } else {
      org <- Mus.musculus::Mus.musculus
    }

    if (is.null(keys)) {
      keys <- get_id(x)()
    }

    out <- rep(NA, length(keys))
    ## check whether is second type id, if is second, repeat by it.
    is.ready <- is_annotation(y)(keys)
    out[is.ready] <- keys[is.ready]

    bool <- !is.ready & is_annotation(x)(keys)
    if (any(bool)) {
      out[bool] <- AnnotationDbi::mapIds(org,
        keys = keys[bool],
        keytype = x,
        column = y,
        multiVals = multiVals
      )
    }
    out
  }
}


#' Search genes annotation in \code{OrganismDb} object
#'
#' Search keytype to more annotation type in \code{OrganismDb} object
#'
#' @param keys Gene id to convert
#' @param species The species of gene. only support Hs or Mm
#' @param multiVals The same with \code{\link[AnnotationDbi]{mapIds}}
#'
#' @seealso
#' \code{\link[AnnotationDbi]{ACCNUM}}
#' \code{\link[Homo.sapiens]{Homo.sapiens}}
#' \code{\link[Mus.musculus]{Mus.musculus}}
#'
#'
#' @name convert_id
#' @rdname convert_id
#' @aliases
#'  alias2symbol alias2entrez alias2ucsc alias2refseq alias2ensembl
#'  symbol2entrez symbol2ensembl symbol2ucsc symbol2refseq symbol2alias
#'  entrez2symbol entrez2ucsc entrez2refseq entrez2alias entrez2entrez
#'  ucsc2symbol ucsc2refseq ucsc2alias ucsc2entrez ucsc2ensembl
#'  refseq2symbol refseq2ucsc refseq2alias refseq2entrez refseq2ensembl
#'  ensembl2symbol ensembl2refseq ensembl2ucsc ensembl2alias ensembl2entrez
#'
#' @return A data frame with annotation info
#'
#' @examples
#' \dontrun{
#'
#' symbol2entrez("TP53", species = "Hs")
#'
#' }
#'
NULL

#' @rdname convert_id
#'
#' @export
#'
alias2symbol <- convert_id("ALIAS", "SYMBOL")

#' @rdname convert_id
#'
#' @export
#'
alias2entrez <- convert_id("ALIAS", "ENTREZID")

#' @rdname convert_id
#'
#' @export
#'
alias2ucsc <- convert_id("ALIAS", "UCSCKG")

#' @rdname convert_id
#'
#' @export
#'
alias2refseq <- convert_id("ALIAS", "REFSEQ")

#' @rdname convert_id
#'
#' @export
#'
alias2ensembl <- convert_id("ALIAS", "ENSEMBL")

#' @rdname convert_id
#'
#' @export
#'
symbol2entrez <- convert_id("SYMBOL", "ENTREZID")

#' @rdname convert_id
#'
#' @export
#'
symbol2ensembl <- convert_id("SYMBOL", "ENSEMBL")

#' @rdname convert_id
#'
#' @export
#'
symbol2refseq <- convert_id("SYMBOL", "REFSEQ")

#' @rdname convert_id
#'
#' @export
#'
symbol2ucsc <- convert_id("SYMBOL", "UCSCKG")

#' @rdname convert_id
#'
#' @export
#'
symbol2alias <- convert_id("SYMBOL", "ALIAS")

#' @rdname convert_id
#'
#' @export
#'
entrez2symbol <- convert_id("ENTREZID", "SYMBOL")

#' @rdname convert_id
#'
#' @export
#'
entrez2ensembl <- convert_id("ENTREZID", "ENSEMBL")

#' @rdname convert_id
#'
#' @export
#'
entrez2refseq <- convert_id("ENTREZID", "REFSEQ")

#' @rdname convert_id
#'
#' @export
#'
entrez2ucsc <- convert_id("ENTREZID", "UCSCKG")

#' @rdname convert_id
#'
#' @export
#'
entrez2alias <- convert_id("ENTREZID", "ALIAS")

#' @rdname convert_id
#'
#' @export
#'
refseq2entrez <- convert_id("REFSEQ", "ENTREZID")

#' @rdname convert_id
#'
#' @export
#'
refseq2symbol <- convert_id("REFSEQ", "SYMBOL")

#' @rdname convert_id
#'
#' @export
#'
refseq2ensembl <- convert_id("REFSEQ", "ENSEMBL")

#' @rdname convert_id
#'
#' @export
#'
refseq2ucsc <- convert_id("REFSEQ", "UCSCKG")

#' @rdname convert_id
#'
#' @export
#'
refseq2alias <- convert_id("REFSEQ", "ALIAS")

#' @rdname convert_id
#'
#' @export
#'
ucsc2symbol <- convert_id("UCSCKG", "SYMBOL")

#' @rdname convert_id
#'
#' @export
#'
ucsc2ensembl <- convert_id("UCSCKG", "ENSEMBL")

#' @rdname convert_id
#'
#' @export
#'
ucsc2refseq <- convert_id("UCSCKG", "REFSEQ")

#' @rdname convert_id
#'
#' @export
#'
ucsc2entrez <- convert_id("UCSCKG", "ENTREZ")

#' @rdname convert_id
#'
#' @export
#'
ucsc2alias <- convert_id("UCSCKG", "ALIAS")
