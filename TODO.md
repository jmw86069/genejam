# TODO for genejam


## Fix bug when colname "intermediate" is supplied as input

### Reproducible example (0.0.13.900 or less)

```
idf <- data.frame(Gene=c("MINA", "MINA", "GABRR3", "GABRR3"),
   intermediate=c("", "84864", "", "200959"))
idf
genejam::freshenGenes2(idf)
```

### Output

```
    Gene intermediate intermediate_v1 SYMBOL                                                               GENENAME
1   MINA                        84864                                                                              
2   MINA        84864           84864  RIOX2                                                  ribosomal oxygenase 2
3 GABRR3                       200959                                                                              
4 GABRR3       200959          200959 GABRR3 gamma-aminobutyric acid type A receptor rho3 subunit (gene/pseudogene)
```

It creates a new column `"intermediate_v1"` then does not use it.

### Desired outcome

Values should be left as-is in `"intermediate"` and new values
should be used to fill in holes.

Ideally, user can specify the `intermediate_colname` so its
values are recognized, potentially skipping the first step
in `freshenGenes()` if there are no remaining colnames
in the input data.

An argument `intermediate_colname` would make the process
a bit more transparent.

Also, the colname `intermediate_colname` should not be split
into multiple columns as is done with input data when it
contains delimited values.

### New workflow (version 0.0.14.900)

```
idf <- data.frame(Gene=c("MINA", "MINA", "GABRR3", "GABRR3", "APOE"),
   ENTREZID=c("", "84864", "", "200959", "84864"))
idf
genejam::freshenGenes2(idf, intermediate="ENTREZID")
genejam::freshenGenes2(idf, intermediate="ENTREZID", handle_multiple="first_hit")
genejam::freshenGenes2(idf, intermediate="ENTREZID", handle_multiple="first_try")
genejam::freshenGenes2(idf, intermediate="ENTREZID", handle_multiple="all")
genejam::freshenGenes2(idf, intermediate="ENTREZID", handle_multiple="best_each")
genejam::freshenGenes2(idf, intermediate="ENTREZID", include_source=TRUE)
```

Trickier example showing mix of input types:

```
idf <- data.frame(Gene=c("MINA", "RIOX2", "GABRR3", "GABRR3", "RIOX2,GABRR3", ""),
   ENTREZID=c("", "84864", "", "200959", "", "84864,200959"))
idf
genejam::freshenGenes2(idf, intermediate="ENTREZID")
genejam::freshenGenes2(idf, intermediate="ENTREZID", include_source=TRUE)
```

Final example showing input with only the `intermediate` column:

```
idf <- data.frame(Gene=c("MINA", "RIOX2", "GABRR3", "GABRR3", "RIOX2,GABRR3", ""),
   ENTREZID=c("", "84864", "", "200959", "", "84864,200959"))
idf
genejam::freshenGenes2(idf[,"ENTREZID",drop=FALSE], intermediate="ENTREZID")
```



## Allow ENTREZID or EG as optional input column(s)

I bring this issue to the top on 30mar2021.

### Problem space

It is not currently possible to provide ENTREZID values at input,
since these values are expected to be in the `"intermediate"` column
after the first step in the process. There is no annotation that
takes ENTREZID and returns ENTREZID.

A close by-pass is to supply data with colname `"intermediate"`,
which works but breaks the actual use of `"intermediate"`.

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
   

