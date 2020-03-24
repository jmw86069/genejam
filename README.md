
<!-- README.md is generated from README.Rmd. Please edit that file -->

Genejam is intended to freshen gene annotations to the current official
standard. It is particularly useful when comparing genes in two
datasets, for example when those datasets may not be using the same gene
symbols to represent equivalent genes.

## Installation

Install using the `remotes` package:

`remotes::install_github("jmw86069/genejam")`

Note: It is recommended not to use the `devtools` package to install
Github packages, mostly because the `devtools` package has many more
components than are required for installation. Instead, `devtools`
includes all components needed to develop R packages, beyond the scope
of installing one such R package.

## Basic usage

The simplest use is to supply a set of gene symbols:

``` r
genejam::freshenGenes(c("APOE", "APOA", "HIST1H1C"))
```

<div class="kable-table">

| input    | intermediate | SYMBOL |
| :------- | :----------- | :----- |
| APOE     | 348          | APOE   |
| APOA     | 4018         | LPA    |
| HIST1H1C | 3006         | H1-2   |

</div>

For slightly more detail, you can edit the argument `final` to include
things like `"SYMBOL"` (default), `"GENENAME"`, `"ALIAS"`, `"ACCNUM"`,
and more.

``` r
genejam::freshenGenes(c("APOE", "APOA", "HIST1H1C"),
   final=c("SYMBOL", 
      "GENENAME",
      "ALIAS"))
```

<div class="kable-table">

| input    | intermediate | SYMBOL | GENENAME                            | ALIAS                             |
| :------- | :----------- | :----- | :---------------------------------- | :-------------------------------- |
| APOE     | 348          | APOE   | apolipoprotein E                    | AD2,APO-E,ApoE4,APOE,LDLCQ5,LPG   |
| APOA     | 4018         | LPA    | lipoprotein(a)                      | AK38,APOA,LP,LPA                  |
| HIST1H1C | 3006         | H1-2   | H1.2 linker histone, cluster member | H1-2,H1.2,H1C,H1F2,H1s-1,HIST1H1C |

</div>

## Comments on usefulness

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
