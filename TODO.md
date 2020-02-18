# TODO for genejam

## Make searches case-insensitive

* This change requires something like `better_mget()` or
perhaps more properly `imget()`.
* The basic mechanism:

```{r, imget}
input <- c("Apoe", "h1-3", "Hist1H1C");
keys <- ls(org.Hs.egSYMBOL2EG)
keymatch <- match(tolower(input), tolower(keys))
keysfound <- !is.na(keymatch);
values <- as.list(jamba::nameVector(rep(NA, length(input)), input))
valuesfound <- mget(keys[keymatch[keysfound]], org.Hs.egSYMBOL2EG, ifnotfound=NA)
values[keysfound] <- valuesfound;
values;
```

After some testing, it appears the slow step is retrieving all
the keys, which is also performed each time the data is queried,
compounding the problem.

Current recommendation is to use `get_anno_db(..., ignore.case=TRUE)`
which returns an environment where all keys are forced to lowercase.
In this case, the environment consumes real memory, and takes longer
to prepare. However once prepared, the environment can be passed
directory to `freshenGenes()` for rapid operation.

## Allow ENTREZID or EG as optional input column(s)

## Problem space

It is not currently possible to provide ENTREZID values at input,
since these values are expected to be in the `"intermediate"` column
after the first step in the process. There is no annotation that
takes ENTREZID and returns ENTREZID.

* Allow `"ENTREZID"` or `"EG"` as an input column, some mechanism to
allow input data to have the intermediate value, without
having to query an annotation.
* It would essentially mark one or more columns with something like 
`"AsIs"` and the values would be delimited and used to populate
the `"intermediate"` column without query.
* There could be a special-case `ann_lib="ASIS"` that would allow
values in the `"AsIs"` columns to pass through, which would give
a mechanism to try querying another column first, then using the
`"ASIS"` column as a secondary option.

### Suggested mechanism

* `asis_colnames=c("ENTREZID","EG")` would define colnames `c("ENTREZID","EG")`
as "ASIS" which means they are not to be queried.
* `ann_lib=c("ASIS", "SYMBOL2EG", "ALIAS2EG")` would process in order:

    1. take all values from the `asis_colnames` and delimit them
    to populate the `"intermediate"` column.
    2. query `"SYMBOL2EG"` using all columns except `asis_colnames`.
    3. query `"ALIAS2EG"` using all columns except `asis_colnames`.

* `ann_lib=c("SYMBOL2EG", "ASIS", "ALIAS2EG")` would process in order:

   1. query `"SYMBOL2EG"` using all columns except `asis_colnames`.
   2. Populate `"intermediate"` with delimited values from the
   `asis_colnames` using rules defined by `handle_multiple`. That
   is, when `handle_multiple="first_try"` then only rows with
   empty values in `"intermediate" column will be populated.
   3. query `"ALIAS2EG"` using all columns except `asis_colnames`.
   

