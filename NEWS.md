# Version 0.0.16.900

## dependency update

* bumped the version dependency on jamba to `0.0.87.900`, due
to a bug introduced in `jamba::cPaste()` that adversely affected
output for accessions that require multiple rounds of querying.

# Version 0.0.15.900

## updates to existing functions

* `freshenGenes()` new argument `ignore.case` which calls
`genejam::imget()` as a drop-in replacement for `mget()`.
The process was improved by calling `AnnotationDbi::keys()`
instead of `AnnotationDbi::ls()`, and this change is at
least an order of magnitude faster. In brief benchmarks,
using `ignore.case=TRUE` adds roughly 0.1 seconds per
annotation in `try_list`, but otherwise is the same
speed regardless the number of input entries.


# Version 0.0.14.900

## updates to existing functions

* `freshenGenes()` was updated so the input can contain
intermediate values, for example ENTREZID values.
In fact, now input can contain a mixture of gene symbols,
ENTREZID intermediate values, and it will fill in the
holes accordingly.

   * new argument `intermediate` to define the colname that
   contains the intermediate values, most commonly EG which
   are Entrez gene ID.
   * Values are propagated in `intermediate` except when
   `handle_multiple="first_hit"` any existing value in
   `intermediate` is used with no further processing.
   All other `handle_multiple` will combine entries
   into `intermediate`.

## new functions

* `is_empty()` is a small helper function to determine which
entries in a vector are either `NA` or `""`.


# Version 0.0.13.900

## new functions

These two new functions are convenience functions.
I often find myself wanting the gene symbol and
long gene name, so now `freshenGenes2()` does that by default.
To add gene aliases, use `freshenGenes3()`.

* `freshenGenes2()` is a simple extension to `freshenGenes()` that
has `"SYMBOL", "GENENAME"` in the output by default.
* `freshenGenes3()` is a simple extension to `freshenGenes()` that
has `"SYMBOL", "GENENAME", "ALIAS"` in the output by default.

# Version 0.0.12.900

## Changes to existing functions

* `get_anno_db()` logic to check for reciprocal annotation names
was updated to cover more scenarios. Specifically, `"org.Hs.egUNIPROT2EG"`
is properly recognized, it previously was not being recognized by
the reciprocal `"org.Hs.egUNIPROT"` and therefore was being skipped.

# Version 0.0.11.900

## Changes to existing functions

* To prepare for a wider release, I decided to rename (!) some
arguments, to have snake_case instead of camelCase for consistency.
I heard myself complaining about my own package, "Why are some arguments
camelCase and others are snake_case? Pick one!" I complain with a
smile on my face, but still it's a fair point.

    * `finalList` is now `final`
    * `tryList` is now `try_list`
    * `annLib` is now `ann_lib`
    
I suppose I should probably rename `freshenGenes()` to `freshen_genes()`.

## Changes to existing functions

* `get_anno_db()` new argument `ignore.case` which will build an
environment where all keys are converted to lowercase. Ultimately,
this option incurs the lowest performance hit, since the keys
only need to be converted once, then the environment can be used
repeatedly with native `mget()` functions.

## New function

* `imget()` case-insensitive `mget()` -- however once I tested it,
I realized this mechanism is fairly slow when using a fairly large
annotation object. Also, if querying the same data multiple times
using `imget()`, there is no re-use and the cost is incurred each
operation -- very much not ideal. This function will likely be
retired soon.


# Version 0.0.10.900

## New functions

* `better_exists()` and `better_get()` which are (not so humbly)
improved versions of `base::exists()` and `base::get()` respectively.
Their sole benefit is to recognize a package prefix in an object name,
so things like `better_exists("base::get")` will return TRUE since
that object does exist; and subsequently `better_get("base::get")` will
return that object. These functions are mostly useful when using
annotation package prefixes such as `"KEGG.db::KEGGPATHID2NAME"`.
I needed a simple way to test if it exists.
`better_exists()` also allows multiple input values.

## Changes to existing functions

* `get_anno_db()` now calls `better_exists()` and `better_get()` which
allows using a package prefix with annotation names.
* `get_anno_db()` argument `"revmap_suffix"` now allows multiple
possible values, it cycles through each until it finds a match, otherwise
returns `NULL`. Some annotations use suffix `"2ENTREZID"` instead of
`"2EG"`, and still others use `"2NAME"`. I'm sure there will be others.

# Version 0.0.9.900

## Changes to existing functions

* `freshenGenes()` new argument `handle_multiple="best_each"` which
returns the best first try for each delimited entry in each input
row. For example `c("APOE","APOA")` will match `"APOE"` as an
authoritatice gene symbol, but `"APOA"` is matched as an alias
to the new gene symbol `"LPA"`. The output will be `"APOA,LPA"`.
Note that output will contain unique entries delimited, but they
will not be sorted.

# Version 0.0.8.900

The next version will have `handle_multiple="best_each"` which
will find the best match for each entry in a set of delimited
gene symbols. Most useful for something like pathway enrichment
results, where the goal is to retain all possible genes, yet
each gene may require a different type of annotation to
find a match. See `TODO.md` for details.

## Updates to existing functions

* `freshenGenes()` includes a new example showing how to recognize
Affymetrix probesets by using a custom search library.
* `freshenGenes()` handles multiple annotation libraries,
mostly in the form of fully described annotation names,
such as `"org.Hs.egSYMBOL"` and `"hgu133plus2ENTREZID"`.
* `freshenGenes()` new option `empty_rule="na"` which will
replace empty entries with `NA`. Other options `empty_rule="blank"`
replaces with `""`, and `empty_rule="original"` replaces empty
entries in the first output column with the original entry in
the first input column.
* `get_anno_db()` now returns``NULL` when an annotation is
not found, instead of throwing an error. This change allows
the calling function to skip missing annotation gracefully
without using `tryCatch()` to catch the error.

# Version 0.0.7.900

## Updates to existing functions

* `freshenGenes()` now properly ignores NA values without throwing
an exception. NA values are left as-is and returned as NA in the
final output.

# Version 0.0.6.900

## Updates to existing functions

* `freshenGenes()` new argument `protect_inline_sep` helps to
prevent splitting single values that may include the same sep
character, for example not splitting `"H4 clustered histone 10, pseudogene"`
into `"H4 clustered histone 10"` and "`pseudogene"`. Also, the
handling of `finalList` uses `sep` as the `split` since that `sep`
is known to have been used in creating the intermediate values,
therefore it should be consistent in the final step. This subtle
change helps allow a more general split pattern in the first
step, such as `"[, ]+"` which splits at comma and/or space,
without splitting at spaces in subsequent steps.

## New functions

* `lgsub()` is a simple wrapper around `gsub()` except that it
operates on character vectors inside a `list`. This function will
probably be moved into the `jamba` package in the near future.
In fact, there will probably be a small family of list-friendly
functions: lgrep(), lvigrep(), lvgrep(), lstrsplit().

# Version 0.0.5.900

## Updates to existing functions

* `get_anno_db()` was updated to handle reverse-map annotations,
for example requesting `org.Hs.egALIAS` and deriving it from
`org.Hs.egALIAS2EG` using `AnnotationDbi::revmap()`.

# Version 0.0.4.900

## Updates to existing functions

* Another update to resolve cascade `stringsAsFactors` bugs.
Least favorite default option in R.

# Version 0.0.3.900

## Updates to existing functions

* Force input `x` to character in `freshenGenes()`.

# Version 0.0.2.900

Note that genejam requires one Bioconductor annotation package,
usually `org.Hs.eg.db` but can be any valid organism, such
as `org.Mm.eg.db` for mouse, or `org.Rn.eg.db` for rat.

## Bug fixes

* Attempt to fix rare cases of NA values in `mget()`
by using `jamba::rmNA()`.

