
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geneInfo

<!-- badges: start -->
<!-- badges: end -->

The goal of geneInfo is to quick get gene info from Bioconductor
`OrganismDb`

## Installation

You can install the development version of geneInfo like so:

``` r
remotes::install_github('Moonerss/geneInfo')
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(geneInfo)
## basic example code

get_gene_info('TP53')
#> 
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:many mapping between keys and columns
#> 
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#>   SYMBOL ENTREZID          GENENAME       GENETYPE         GO
#> 1   TP53     7157 tumor protein p53 protein-coding GO:0000122
#>                                                      GOTERM TXCHROM
#> 1 negative regulation of transcription by RNA polymerase II   chr17
get_exons_info('TP53')
#> 'select()' returned 1:many mapping between keys and columns
#> # A tibble: 1 × 7
#> # Groups:   exon_chr [1]
#>   exon_chr gene  n_allowed_exons n_total_exons exon_start exon_end gene_length
#>   <chr>    <chr>           <int>         <int>      <int>    <int>       <dbl>
#> 1 chr17    TP53               12            21    7565097  7590868       25772
allowedExons('TP53')
#> 'select()' returned 1:many mapping between keys and columns
#> TP53 
#>   12
totalExons('TP53')
#> 'select()' returned 1:many mapping between keys and columns
#> TP53 
#>   21
geneLength('TP53')
#> 'select()' returned 1:many mapping between keys and columns
#>  TP53 
#> 25772
is_symbol('TP53', species = 'Hs')
#> [1] TRUE
is_symbol(c('TP53', 'SJDSJ'))
#> [1]  TRUE FALSE
get_geneid_type(c('TP53', 'ADSL'), species = 'Hs')
#>     TP53     ADSL 
#> "SYMBOL" "SYMBOL"
head(get_ucsc())
#> [1] "ENST00000596924.1" "ENST00000598345.1" "ENST00000600966.1"
#> [4] "ENST00000595014.1" "ENST00000263100.8" "ENST00000596636.1"
search_symbol('TP53', species = 'Hs')
#> 'select()' returned 1:1 mapping between keys and columns
#>   SYMBOL ENTREZID TXCHROM
#> 1   TP53     7157   chr17
symbol2entrez("TP53", species = "Hs")
#> 'select()' returned 1:1 mapping between keys and columns
#> [1] "7157"
```
