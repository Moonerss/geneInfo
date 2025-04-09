# Build function to get info
#
# Build function based on some argument
#
# @param keytype keytype in \code{OrganismDb} object
#
# @return return a function object to get specific keytype
#
#
get_id <- function(keytype) {
  .get_id <- function(species = c('Hs', 'Mm')) {
    species <- match.arg(species)
    if (species == 'Hs') {
      org <- Homo.sapiens::Homo.sapiens
    } else {
      org <- Mus.musculus::Mus.musculus
    }
    AnnotationDbi::keys(org, keytype=keytype)
  }
}


#'
#' Get keytype from \code{OrganismDb} object
#'
#' Get different keytype info from \code{OrganismDb} object
#'
#' @param species The species of gene. only support Hs or Mm
#'
#' @name get_id
#' @rdname get_id
#'
#' @seealso
#' \code{\link[AnnotationDbi]{ACCNUM}}
#' \code{\link[Homo.sapiens]{Homo.sapiens}}
#' \code{\link[Mus.musculus]{Mus.musculus}}
#'
#' @examples
#' \dontrun{
#'
#' get_ucsc(species = 'Hs')
#'
#' get_refseq(species = 'Hs')
#'
#' get_symbol(species = 'Hs')
#'
#' get_ensembl(species = 'Hs')
#'
#' get_entrez(species = 'Hs')
#'
#' get_genebank(species = 'Hs')
#' }
#'
NULL


#'
#' @rdname get_id
#'
#' @export
#'
get_ucsc <- get_id(keytype = 'UCSCKG')


#' @rdname get_id
#'
#' @export
#'
get_refseq <- get_id(keytype = 'REFSEQ')

#' @rdname get_id
#'
#' @export
#'
get_symbol <- get_id(keytype = 'SYMBOL')

#' @rdname get_id
#'
#' @export
#'
get_ensembl <- get_id(keytype = 'ENSEMBL')

#' @rdname get_id
#'
#' @export
#'
get_entrez <- get_id(keytype = 'ENTREZID')

#' @rdname get_id
#'
#' @export
#'
get_genebank <- get_id(keytype = 'ACCNUM')


# Search genes annotation
#
# Search keytype to more annotation type
#
# @param keytype the annotation keytype
#
# @seealso
# \code{\link[AnnotationDbi]{ACCNUM}}
# \code{\link[Homo.sapiens]{Homo.sapiens}}
# \code{\link[Mus.musculus]{Mus.musculus}}
#
search_genes <- function(keytype) {
  .search_genes <- function(keys=NULL, species = c('Hs', 'Mm'), columns=NULL, add.columns=NULL, remove.columns=NULL) {
    species <- match.arg(species)
    if (species == 'Hs') {
      org <- Homo.sapiens::Homo.sapiens
    } else {
      org <- Mus.musculus::Mus.musculus
    }
    if (is.null(columns)) {
      columns <- c('SYMBOL','ENTREZID','TXCHROM')
    }
    columns <- setdiff(columns,keytype)
    columns <- setdiff(columns,remove.columns)
    columns <- union(columns,add.columns)

    if (is.null(keys)) {
      keys <- AnnotationDbi::keys(org, keytype=keytype)
    }

    AnnotationDbi::select(org,keys=keys,keytype=keytype,columns=columns)
  }
}

#' Search genes annotation in \code{OrganismDb} object
#'
#' Search keytype to more annotation type in \code{OrganismDb} object
#'
#' @param keys Gene id to search
#' @param species The species of gene. only support Hs or Mm
#' @param columns The column to return. If NULL, return three column: SYMBOL, ENTREZID, and TXCHROM
#' @param add.columns The column to add in the result, the column must in the \code{OrganismDb} object
#' @param remove.columns The column to remove in the result, the column must in the result.
#'
#' @seealso
#' \code{\link[AnnotationDbi]{ACCNUM}}
#' \code{\link[Homo.sapiens]{Homo.sapiens}}
#' \code{\link[Mus.musculus]{Mus.musculus}}
#'
#'
#' @name search_genes
#' @rdname search_genes
#' @aliases search_entrez search_symbol search_ensembl search_refseq search_alias
#'
#' @return A data frame with annotation info
#'
#' @examples
#' \dontrun{
#'
#' search_symbol('TP53', species = 'Hs')
#'
#' }
#'
#'
NULL

#' @rdname search_genes
#'
#' @export
#'
search_symbol <- search_genes('SYMBOL')

#' @rdname search_genes
#'
#' @export
#'
search_ensembl <- search_genes('ENSEMBL')

#' @rdname search_genes
#'
#' @export
#'
search_entrez <- search_genes('ENTREZID')

#' @rdname search_genes
#'
#' @export
#'
search_refseq <- search_genes('REFSEQ')

#' @rdname search_genes
#'
#' @export
#'
search_alias <- search_genes('ALIAS')


