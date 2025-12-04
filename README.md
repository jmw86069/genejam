
Genejam is intended to freshen gene annotations to the current official
standard. It is particularly useful when comparing genes in two
datasets, for example when those datasets may not be using the same gene
symbols to represent equivalent genes.

The notable behavior is that it uses the “best available” source of
annotation cross-reference, with ordered logic.

1.  If ENTREZID is a perfect match, use that. (NCBI Entrez gene ID)
2.  If SYMBOL is a perfect match, use that. (NCBI Entrez gene symbol)
3.  If ACCNUM is a match, it uses that. (Any accession number)
4.  If ALIAS is a match, it uses that. (Other gene aliases, xrefs)

The ‘ALIAS’ option may find multiple genes. Sometimes this is the best
available option, so be it. However, it is only used after trying
ENTREZID, then official gene symbol, then accession numbers first.

This step may be customized, for example it is appropriate to include
other accession number sources such as: UNIPROT, ENSEMBL, ENSEMBLTRANS,
and platform identifiers such as Affymetrix or Illumina PROBE_SET_ID or
PROBEID.

Other important features:

- For any input which matches more than one entry, multiple values are
  returned.
- The query can use case-insensitive steps, which is helpful for mouse
  gene symbols which are mixed-case “Ccl8”. However, there are also
  genes in mouse ‘Gs’ is the gene “greasy” ENTREZID 109550, and
  lowercase ‘gs’ is the gene “glabrous”, ENTREZID 14835. Sometimes the
  specific capitalization is useful, so the default method uses
  `ignore.case=FALSE`.

## Installation

Install using the `remotes` package:

``` r
remotes::install_github("jmw86069/genejam")
```

or with the `pak` package:

``` r
pak::pkg_install("jmw86069/genejam")
```

## freshenGenes()

The simplest use is to supply a set of gene symbols:

``` r
genejam::freshenGenes(c("APOE", "APOA", "HIST1H1C"))
```

<div class="kable-table">

| input    | ENTREZID | SYMBOL |
|:---------|:---------|:-------|
| APOE     | 348      | APOE   |
| APOA     | 4018     | LPA    |
| HIST1H1C | 3006     | H1-2   |

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

| input | ENTREZID | SYMBOL | GENENAME | ALIAS |
|:---|:---|:---|:---|:---|
| APOE | 348 | APOE | apolipoprotein E | AD2,APO-E,ApoE4,APOE,LDLCQ5,LPG |
| APOA | 4018 | LPA | lipoprotein(a) | AK38,APOA,LP,LPA |
| HIST1H1C | 3006 | H1-2 | H1.2 linker histone, cluster member | H1-2,H1.2,H1C,H1F2,H1s-1,HIST1H1C |

</div>

## Convenience function freshenGenes2()

I frequently find myself wanting gene symbol, and the long gene name, so
I created a simple function `freshenGenes2()` that uses default argument
`final=c("SYMBOL", "GENENAME")`:

``` r
genejam::freshenGenes2(c("APOE", "APOA", "HIST1H1C"))
```

<div class="kable-table">

| input    | ENTREZID | SYMBOL | GENENAME                            |
|:---------|:---------|:-------|:------------------------------------|
| APOE     | 348      | APOE   | apolipoprotein E                    |
| APOA     | 4018     | LPA    | lipoprotein(a)                      |
| HIST1H1C | 3006     | H1-2   | H1.2 linker histone, cluster member |

</div>

## Convenience function freshenGenes3()

The other common use case is to include other gene aliases, with the
function `freshenGenes3()`:

``` r
genejam::freshenGenes3(c("APOE", "APOA", "HIST1H1C"))
```

<div class="kable-table">

| input | ENTREZID | SYMBOL | GENENAME | ALIAS |
|:---|:---|:---|:---|:---|
| APOE | 348 | APOE | apolipoprotein E | AD2,APO-E,ApoE4,APOE,LDLCQ5,LPG |
| APOA | 4018 | LPA | lipoprotein(a) | AK38,APOA,LP,LPA |
| HIST1H1C | 3006 | H1-2 | H1.2 linker histone, cluster member | H1-2,H1.2,H1C,H1F2,H1s-1,HIST1H1C |

</div>

## Slightly more advanced

What if you already have Entrez gene ID, and want associated
annotations? The function `freshenGenes()` runs two steps:

1.  convert input to `intermediate` (which is the Entrez gene ID)
2.  convert `intermediate` to the output defined by `final`, for example
    `final=c("SYMBOL")` would create a column `"SYMBOL"`.

In this example, the Entrez gene ID values are in a column `"ENTREZID"`,
so we will use argument `intermediate="ENTREZID"`. In this case, you
already have `"intermediate"`, so you invoke the function with a
`data.frame` with values in a column named `"intermediate"`, and set
`try_list=NULL`.

``` r
df <- data.frame(ENTREZID=c("348", "4018", "3006", "100"));
genejam::freshenGenes2(x=df, intermediate="ENTREZID")
```

<div class="kable-table">

| ENTREZID | SYMBOL | GENENAME                            |
|:---------|:-------|:------------------------------------|
| 348      | APOE   | apolipoprotein E                    |
| 4018     | LPA    | lipoprotein(a)                      |
| 3006     | H1-2   | H1.2 linker histone, cluster member |
| 100      | ADA    | adenosine deaminase                 |

</div>

Similarly, you can provide input with a mixture of gene symbols and
Entrez gene ID values. Shown below is mixed input.

``` r
idf <- data.frame(Gene=c("MINA", "", "GABRR3", "GABRR3", ""),
   ENTREZID=c("", "84864", "", "200959", "200959"))
idf
```

<div class="kable-table">

| Gene   | ENTREZID |
|:-------|:---------|
| MINA   |          |
|        | 84864    |
| GABRR3 |          |
| GABRR3 | 200959   |
|        | 200959   |

</div>

You only need to specify `intermediate="ENTREZID"` as before.

``` r
genejam::freshenGenes2(x=idf, intermediate="ENTREZID")
```

<div class="kable-table">

| Gene   | ENTREZID | SYMBOL | GENENAME                                             |
|:-------|:---------|:-------|:-----------------------------------------------------|
| MINA   | 84864    | RIOX2  | ribosomal oxygenase 2                                |
|        | 84864    | RIOX2  | ribosomal oxygenase 2                                |
| GABRR3 | 200959   | GABRR3 | gamma-aminobutyric acid type A receptor subunit rho3 |
| GABRR3 | 200959   | GABRR3 | gamma-aminobutyric acid type A receptor subunit rho3 |
|        | 200959   | GABRR3 | gamma-aminobutyric acid type A receptor subunit rho3 |

</div>

Notice the values in `"ENTREZID"` are updated based upon the first step
resolution of `"Gene"` values to `"ENTREZID"`. The `"SYMBOL"` and
`"GENENAME"` columns are populated using values in `"ENTREZID"`.

## When is genejam useful?

The official gene nomenclature is updated multiple times per year, which
means one Entrez gene ID may have a different official gene symbol
before and after the update. When comparing data from two experiments,
it is important to use the same gene nomenclature. Otherwise, there will
be differences in results only because the names of some genes are
different.

Most microarray platforms provide gene annotations, which are updated
much less frequently than the official genes. For example, Affymetrix
array “Clariom D Human” was last updated between 2016 and 2018 (this
document was written in 2021.) In order to compare microarray results to
those from literature, biological pathways, or other experiments, the
gene nomenclature needs to be updated to the most current version.

### Rationale for the workflow

In rare cases an official gene symbol is “moved” from one Entrez gene ID
to another, usually when the original Entrez gene ID is deleted. In
these cases, the most reasonable link between an experimental asay and
the targeted gene is the gene symbol. An alternative is to use a
sequence accession number used to design the assay.

As a result, a “best possible” gene annotation strategy is used.

- When Entrez gene ID is available, use it.
- When Entrez official gene symbol is available, use that to determine
  the Entrez gene ID.
- When a sequence accession used for assay design is available, use that
  to determine the associated Entrez gene ID value or values.
- When a gene symbol alias is available, use that to determine the
  Entrez gene ID value or values associated with this alias.

Sometimes an assay measures two genes. The steps in `genejam` are
designed to maintain multi-gene associations where necessary. If one
gene symbol alias is associated with two (or more) genes, then all those
genes will become associated. Note this only happens if a higher
priority association was not already found.

Ultimately the workflow is what I and others have been doing all along,
to assemble the best available gene annotations for a given dataset,
while also leaving behind the fewest possible un-annotated entries. When
one source of annotation fails, try another on missing entries; and so
on.

## Optimizations

The steps used in `genejam` are designed for speed, to the extent that
providing 100,000 rows should return results within a few seconds at
most. Only unique values are queried, and only missing entries are
updated. When multiple values are combined by a delimiter, a highly
optimized method is used.

Lastly, these operations also use optimized mixed-alphanumeric sorting
even in the context of a `list`, so things like `"chr2"` will appear
before `"chr10"` in sort order. Incidentally, a sort step is necessary
so you can compare whether two entries are associated with the same
genes. If one is stored as `"APOE,APOE4"` and the other as
`"APOE4,APOE"`, this comparison fails.

All that to say, I use these functions a lot so I need them to be
reliable and fast. It takes a few seconds just to load the associated
SQLite gene annotation data, and the process used by `freshenGenes()` is
usually substantially faster than that step alone.

## Package Reference

A full online function reference is available via the pkgdown
documentation:

Full genejam command reference: <https://jmw86069.github.io/genejam>
