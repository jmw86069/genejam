# TODO for genejam


## New option `handle_multiple="each_best"`

This option will essentially mimic `handle_multiple="first_try"`
except that when the input data has delimited entries on each line,
each delimited entry will be split and handled individually using
the `"first_try"` logic.

For example, consider an entry with:

> `"APOE,HIST1H2A"`

The try list is `tryList=c("SYMBOL2EG","ACCNUM2EG","ALIAS2EG")` which means
it will try to match:

* the official gene symbol `"SYMBOL2EG"` and stop if it finds a match
* then accession numbers `"ACCNUM2EG"` and stop if it finds a match
* then aliases `"ALIAS2EG"` and stop if it finds a match.

The best entry for `"APOE"` is from `"SYMBOL2EG"`.
The best entry for `"HIST1H2A"` is from `"ALIAS2EG"`.

If we use `handle_multiple="first_try"` with `"APOE,HIST1H2A"` the
output will be `"APOE"` -- since we stop trying after the first
`tryLib` entry is successful.

Instead, we want to split `"APOE,HIST1H2A"` into `"APOE"` and
`"HIST1H2A"`, then match each separately, then combine the results.

When each entry contains aliases of the same gene, use `handle_multiple="first_try"`.
When each entry contains different genes, and you want the best
current annotation for each gene, use `handle_multiple="each_best"`.

### Basic workflow

The general guidance is to use only one input column at a time.

```
require(org.Hs.eg.db);
## Define an interesting input list
input <- c("APOE,HIST1H1C,HIST1H1D,,HIST1H2A", "HIST1H2A","APOA", "APOA,LPA,HIST1H1C");

## Split the input by delimiter
taller_list <- strsplit(input, ",");
## Make a vector and associated factor to split back into a list
taller_vector <- unlist(taller_list);
taller_factor <- rep(factor(seq_along(taller_list)),
   lengths(taller_list));

## run freshenGenes()
taller_freshened <- freshenGenes(taller_vector,
   empty_rule="na",
   handle_multiple="first_try");
taller_freshened;

## optionally review any changes
changed_rows <- subset(unique(taller_freshened), jamba::rmNA(naValue="",input) != jamba::rmNA(naValue="",SYMBOL));
changed_rows;

## Review changes entries in more detail
jamba::mixedSortDF(byCols="SYMBOL", 
   subset(unique(taller_freshened), 
      SYMBOL %in% jamba::rmNA(changed_rows$SYMBOL)))

## Split back into the original vector
taller_freshened_split <- split(taller_freshened$SYMBOL, taller_factor);

## comma-delimit using only unique entries
jamba::cPasteU(taller_freshened_split, na.rm=TRUE);
```

