
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

## freshenGenes()

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

## Convenience function freshenGenes2()

I frequently find myself wanting gene symbol, and the long gene name, so
I created a simple function `freshenGenes2()` that uses default argument
`final=c("SYMBOL", "GENENAME")`:

``` r
genejam::freshenGenes2(c("APOE", "APOA", "HIST1H1C"))
```

<div class="kable-table">

| input    | intermediate | SYMBOL | GENENAME                            |
| :------- | :----------- | :----- | :---------------------------------- |
| APOE     | 348          | APOE   | apolipoprotein E                    |
| APOA     | 4018         | LPA    | lipoprotein(a)                      |
| HIST1H1C | 3006         | H1-2   | H1.2 linker histone, cluster member |

</div>

## Convenience function freshenGenes3()

The other common use case is to include other gene aliases, with the
function `freshenGenes3()`:

``` r
genejam::freshenGenes3(c("APOE", "APOA", "HIST1H1C"))
```

<div class="kable-table">

| input    | intermediate | SYMBOL | GENENAME                            | ALIAS                             |
| :------- | :----------- | :----- | :---------------------------------- | :-------------------------------- |
| APOE     | 348          | APOE   | apolipoprotein E                    | AD2,APO-E,ApoE4,APOE,LDLCQ5,LPG   |
| APOA     | 4018         | LPA    | lipoprotein(a)                      | AK38,APOA,LP,LPA                  |
| HIST1H1C | 3006         | H1-2   | H1.2 linker histone, cluster member | H1-2,H1.2,H1C,H1F2,H1s-1,HIST1H1C |

</div>

## Slightly more advanced

What if you already have Entrez gene ID, and want associated
annotations? The function `freshenGenes()` runs two steps:

1.  convert input to `"intermediate"` (which is the Entrez gene ID)
2.  convert `"intermediate"` to the output defined by `final`, for
    example `final=c("SYMBOL")`.

In this case, you already have `"intermediate"`, so you invoke the
function with a `data.frame` with values in a column named
`"intermediate"`, and set `try_list=NULL`:

``` r
df <- data.frame(intermediate=c("348", "4018", "3006", "100"));
genejam::freshenGenes2(x=df, try_list=NULL)
```

<div class="kable-table">

| intermediate | SYMBOL | GENENAME                            |
| :----------- | :----- | :---------------------------------- |
| 348          | APOE   | apolipoprotein E                    |
| 4018         | LPA    | lipoprotein(a)                      |
| 3006         | H1-2   | H1.2 linker histone, cluster member |
| 100          | ADA    | adenosine deaminase                 |

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
