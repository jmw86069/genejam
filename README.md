
<!-- README.md is generated from README.Rmd. Please edit that file -->

Genejam is intended to freshen gene annotations to the current official
standard. It is particularly useful when comparing genes in two
datasets, for example when those datasets may not be using the same gene
symbols to represent equivalent genes.

``` r
genejam::freshenGenes(c("APOE", "APOA", "HIST1H1C"))
#>      input intermediate SYMBOL
#> 1     APOE          348   APOE
#> 2     APOA         4018    LPA
#> 3 HIST1H1C         3006   H1-2
```

Specifically, in December 2019 the gene nomenclature for the large
histone gene family was changed, so the gene symbol `"HIST1H1C"` from
November 2019 is now called `"H1-2"` in December 2019 and thereafter.
These gene nomenclature updates happen frequently (perhaps monthly), and
apparently any gene symbol can be updated. So when we compare a gene
list created today, to any gene list not created today, we run
`freshenGenes()` on each gene list to ensure they both use the same,
“most current” gene annotation.

## Package Reference

A full online function reference is available via the pkgdown
documentation:

Full genejam command reference: <https://jmw86069.github.io/genejam>
