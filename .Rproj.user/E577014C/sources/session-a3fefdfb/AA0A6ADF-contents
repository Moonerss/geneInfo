#' @title Get Genomic Information
#' @description
#' Get gene info from bioconductor \code{db} package and corresponding \code{Txdb} package
#'
#' @param symbols character vector of HGNC gene symbols
#' @param species The species of gene. only support Hs or Mm
#' @export
#' @rdname get_gene_info
#' @examples
#' \dontrun{
#' get_gene_info('TP53')
#' get_gene_info('Pzp', species = 'Mm')
#' }
#'
get_gene_info <- function(symbols, species = c('Hs', 'Mm')) {
  # see https://davetang.org/muse/2013/12/16/bioconductor-annotation-packages/
  species <- match.arg(species)
  columns = c('ENTREZID', 'GENENAME', 'GENETYPE')
  if (species == 'Hs') {
    orgdb <- org.Hs.eg.db::org.Hs.eg.db
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  } else {
    orgdb <- org.Mm.eg.db::org.Mm.eg.db
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  }
  info = AnnotationDbi::select(orgdb,
                               keys = symbols,
                               keytype = 'SYMBOL',
                               columns = columns)
  info = info[!duplicated(info$SYMBOL), ]
  info$GO = AnnotationDbi::mapIds(orgdb,
                                  keys = info$SYMBOL,
                                  keytype = 'SYMBOL',
                                  column = 'GO')
  info$GOTERM = AnnotationDbi::mapIds(GO.db::GO.db,
                                      keys = info$GO,
                                      keytype = 'GOID',
                                      column = 'TERM')
  info$TXCHROM = AnnotationDbi::mapIds(txdb,
                                       keys = info$ENTREZID,
                                       keytype = 'GENEID',
                                       column = 'TXCHROM')
  info
}
