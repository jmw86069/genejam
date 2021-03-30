

#' Freshen gene annotations using Bioconductor annotation data
#' 
#' Freshen gene annotations using Bioconductor annotation data
#' 
#' This function takes a vector or `data.frame` of gene symbols,
#' and uses Bioconductor annotation methods to find the most current
#' official gene symbol.
#' 
#' The annotation process runs in two basic steps:
#'
#' 1. **Convert the input gene to Entrez gene ID**.
#' 2. **Convert Entrez gene ID to official gene symbol**.
#'
#' ## Step 1. Convert to Entrez gene ID
#'
#' The first step uses an ordered list of annotations,
#' with the assumption that the first match is usually the best,
#' and most specific. By default, the order is:
#'
#' * `"org.Hs.egSYMBOL2EG"` -- almost always 1-to-1 match
#' * `"org.Hs.egACCNUM2EG"` -- mostly a 1-to-1 match
#' * `"org.Hs.egALIAS2EG"` -- sometimes a 1-to-1 match, sometimes 1-to-many
#'
#' When multiple Entrez gene ID values are matched, they are all
#' retained. See argument `handle_multiple` for custom options.
#'
#' ## Step 2. Use Entrez gene ID to return official annotation
#'
#' The second step converts the Entrez gene ID (or multiple IDs)
#' to the official gene symbol, by default using `"org.Hs.egSYMBOL"`.
#'
#' The second step may optionally include multiple annotation types,
#' each of which will be returned. Some common examples:
#'
#' * `"org.Hs.egSYMBOL"` -- official Entrez gene symbol
#' * `"org.Hs.egALIAS"` -- set of recognized aliases for an Entrez gene.
#' * `"org.Hs.egGENENAME"` -- official Entrez long gene name
#'
#' For each step, the annotation matched can be returned, as an audit
#' trail to see which annotation was available for each input entry.
#' 
#' Note that if the input data already contains Entrez gene ID
#' values, you can define that colname with argument `intermediate`.
#' 
#' @return `data.frame` with one or more columns indicating the input
#' data, then a column `"intermediate"` containing the Entrez gene ID
#' that was matched, then one column for each item in `final`,
#' by default `"SYMBOL"`.
#' 
#' @family genejam
#' 
#' @param x character vector or `data.frame` with one or most columns
#'    containing gene symbols.
#' @param ann_lib character vector indicating the name or names of the
#'    Bioconductor annotation library to use when looking up
#'    gene nomenclature.
#' @param try_list character vector indicating one or more names of
#'    annotations to use for the input gene symbols in `x`. The
#'    annotation should typically return the Entrez gene ID, usually
#'    given by `'2EG'` at the end of the name. For example `SYMBOL2EG`
#'    will be used with ann_lib `"org.Hs.eg.db"` to produce annotation
#'    name `"org.Hs.egSYMBOL2EG"`. Note that when the `'2EG'` form of
#'    annotation does not exist (or another suitable suffix defined in 
#'    argument `"revmap_suffix"` in `get_anno_db()`), it will be derived
#'    using `AnnotationDbi::revmap()`. For example if `"org.Hs.egALIAS"`
#'    is requested, but only `"org.Hs.egALIAS2EG"` is available, then
#'    `AnnotationDbi::revmap(org.Hs.egALIAS2EG)` is used to create the
#'    equivalent of `"org.Hs.egALIAS"`.
#' @param final character vector to use for the final conversion
#'    step. When `final` is `NULL` no conversion is performed.
#'    When `final` contains multiple values, each value is returned
#'    in the output. For example, `final=c("SYMBOL","GENENAME")` will
#'    return a column `"SYMBOL"` and a column `"GENENAME"`.
#' @param split character value used to separate delimited values in `x`
#'    by the function `base::strsplit()`. The default will split values
#'    separated by comma `,` semicolon `;` or forward slash `/`, and will
#'    trim whitespace before and after these delimiters.
#' @param sep character value used to concatenate multiple entries in
#'    the same field. The default `sep=","` will comma-delimit multiple
#'    entries in the same field.
#' @param handle_multiple character value indicating how to handle multiple
#'    values: `"first_hit"` will query each column of `x` until it finds the
#'    first possible returning match, and will ignore all subsequent possible
#'    matches for that row in `x`. For example, if one row in `x` contains
#'    multiple values, only the first match will be used. `"first_try"`
#'    will return the first match from `try_list` for all columns in `x`
#'    that contain a match. For example, if one row in `x` contains two
#'    values, the first match from `try_list` using one or both columns in
#'    `x` will be maintained. Subsequent entries in `try_list` will not be
#'    attempted for rows that already have a match. `"all"` will return all
#'    possible matches for all entries in `x` using all items in `try_list`.
#' @param empty_rule character value indicating how to handle entries which
#'    did not have a match, and are therefore empty: `"original"` will use
#'    the original entry as the output field; `"empty"` will leave the
#'    entry blank.
#' @param include_source logical indicating whether to include a column
#'    that shows the colname and source matched. For example, if column
#'    `"original_gene"` matched `"SYMBOL2EG"` in `"org.Hs.eg.db"` there
#'    will be a column `"found_source"` with value
#'    `"original_gene.org.Hs.egSYMBOL2EG"`.
#' @param protect_inline_sep logical indicating whether to
#'    protect inline characters in `sep`, to prevent them from
#'    being used to split single values into multiple values.
#'    For example, `"GENENAME"` returns the full gene name, which
#'    often contains comma `","` characters. These commas do
#'    not separate multiple separate values, so they should not be
#'    used to split a string like `"H4 clustered histone 10, pseudogene"`
#'    into two strings `"H4 clustered histone 10"` and `"pseudogene"`.
#' @param intermediate `character` string with colname in `x` that
#'    contains intermediate values. These values are expected from output
#'    of the first step in the workflow, for example `"SYMBOL2EG"`
#'    returns Entrez gene values, so if the input `x` already contains
#'    some of these values in a column, assign that colname to
#'    `intermediate`.
#' @param verbose logical indicating whether to print verbose output.
#' 
#' @examples
#' if (suppressPackageStartupMessages(require(org.Hs.eg.db))) {
#'    cat("\nBasic usage\n");
#'    print(freshenGenes(c("APOE", "CCN2", "CTGF")));
#' }
#' 
#' if (suppressPackageStartupMessages(require(org.Hs.eg.db))) {
#'    ## Optionally show the annotation source matched
#'    cat("\nOptionally show the annotation source matched\n");
#'    print(freshenGenes(c("APOE", "CCN2", "CTGF"), include_source=TRUE));
#' }
#' 
#' if (suppressPackageStartupMessages(require(org.Hs.eg.db))) {
#'    ## Show comma-delimited genes
#'    cat("\nInput genes are comma-delimited\n");
#'    print(freshenGenes(c("APOE", "CCN2", "CTGF", "CCN2,CTGF")));
#' }
#' 
#' if (suppressPackageStartupMessages(require(org.Hs.eg.db))) {
#'    ## Optionally include more than SYMBOL in the output
#'    cat("\nCustom output to include SYMBOL, ALIAS, GENENAME\n");
#'    print(freshenGenes(c("APOE", "HIST1H1C"),
#'       final=c("SYMBOL", "ALIAS", "GENENAME")));
#' }
#' 
#' if (suppressPackageStartupMessages(require(org.Hs.eg.db))) {
#'    ## More advanced, match affymetrix probesets
#'    if (suppressPackageStartupMessages(require(hgu133plus2.db))) {
#'       cat("\nAdvanced example including Affymetrix probesets.\n");
#'       print(freshenGenes(c("227047_x_at","APOE","HIST1H1D","NM_003166,U08032"),
#'          include_source=TRUE,
#'          try_list=c("hgu133plus2ENTREZID","REFSEQ2EG","SYMBOL2EG","ACCNUM2EG","ALIAS2EG"),
#'          final=c("SYMBOL","GENENAME")))
#'    }
#' }
#' 
#' @export
freshenGenes <- function
(x,
 ann_lib=c("","org.Hs.eg.db"),
 try_list=c("SYMBOL2EG", "ACCNUM2EG", "ALIAS2EG"),
 final=c("SYMBOL"),
 split="[ ]*[,/;]+[ ]*",
 sep=",",
 handle_multiple=c("first_try", "first_hit", "all", "best_each"),
 empty_rule=c("empty", "original", "na"),
 include_source=FALSE,
 protect_inline_sep=TRUE,
 intermediate="intermediate",
 verbose=FALSE,
 ...)
{
   ###
   handle_multiple <- match.arg(handle_multiple);
   empty_rule <- match.arg(empty_rule);
   if (length(ann_lib) == 0) {
      ann_lib <- "";
   }
   test_ann_lib <- setdiff(ann_lib, c(NA,""));
   if (length(test_ann_lib) > 0) {
      for (test_lib in test_ann_lib) {
         if (!suppressPackageStartupMessages(require(test_lib, character.only=TRUE))) {
            stop(paste0("Package '", test_lib, "' is not available."));
         }
      }
   }
   if (is.atomic(x)) {
      x <- data.frame(input=as.character(x),
         stringsAsFactors=FALSE,
         check.names=FALSE);
   }
   if (length(colnames(x)) == 0) {
      colnames(x) <- jamba::makeNames(rep("input", ncol(x)));
   }
   
   ## colnames_x are the colnames(x) that are not intermediate
   intermediate_source <- paste0(intermediate, "_source");
   colnames_x <- setdiff(colnames(x),
      c(intermediate, intermediate_source));
   ## ncol_x is the number of columns that are not intermediate
   ncol_x <- length(colnames_x);
   
   ## handle_multiple="best_each"
   if ("best_each" %in% handle_multiple) {
      if (verbose) {
         jamba::printDebug("freshenGenes(): ",
            "handle_multiple:",
            handle_multiple);
      }
      if (ncol_x == 1) {
         ## Check for delimited values
         if (length(split) > 0 && nchar(split) > 0 && jamba::igrepHas(split, x[[colnames_x]])) {
            ## Split the input by delimiter
            taller_list <- strsplit(x[[colnames_x]], split);
            # 30mar2021: fill NULL with "" so the original empty entry is not lost
            taller_list[lengths(taller_list) == 0] <- "";
            ## Make a vector and associated factor to split back into a list
            taller_idx <- rep(seq_along(taller_list),
               lengths(taller_list));
            taller_factor <- rep(factor(seq_along(taller_list)),
               lengths(taller_list));
            taller_x <- x[taller_idx,,drop=FALSE];
            taller_x[[colnames_x]] <- unlist(taller_list);
            #taller_vector <- unlist(taller_list);
         } else {
            taller_x <- x;
            #taller_vector <- x[[colnames_x]];
            taller_factor <- NULL;
         }
         ## run freshenGenes()
         ## empty_rule="na" here so blank entries will get dropped
         ## then we can replace as needed later
         taller_freshened <- freshenGenes(x=taller_x,
            ann_lib=ann_lib,
            try_list=try_list,
            final=final,
            split=split,
            sep=sep,
            handle_multiple="first_try",
            empty_rule="na",
            include_source=include_source,
            protect_inline_sep=protect_inline_sep,
            intermediate=intermediate,
            verbose=FALSE);
         ## Split back into the original vector
         if (length(taller_factor) == 0) {
            return(taller_freshened);
         }
         # 30mar2021: this chunk was ignored, commenting out
         #final_colnames <- jamba::provigrep(c(final,
         #   #"^intermediate",
         #   intermediate,
         #   "_source$"),
         #   colnames(taller_freshened));
         final_colnames <- colnames(taller_freshened);
         x_new <- do.call(cbind, lapply(jamba::nameVector(final_colnames), function(i){
            taller_freshened_split <- split(taller_freshened[[i]],
               taller_factor);
            ## comma-delimit using only unique entries
            idf <- data.frame(check.names=FALSE,
               stringsAsFactors=FALSE,
               output=jamba::cPasteU(taller_freshened_split,
                  sep=sep,
                  na.rm=TRUE));
            colnames(idf) <- i;
            idf;
         }));
         return(x_new);
      }
   }
   
   ## Expand columns containing delimited values if necessary
   # This step makes multiple values appear in separate columns
   # on the same row.
   if (length(split) > 0 && nchar(split) > 0) {
      x <- data.frame(stringsAsFactors=FALSE,
         check.names=FALSE,
         do.call(cbind,
            lapply(jamba::nameVector(colnames(x)), function(i){
               ix <- as.character(x[[i]]);
               # only split delimited values when it is not intermediate
               # note that empty entries are filled with ""
               if (i %in% colnames_x) {
                  if (jamba::igrepHas(split, ix)) {
                     ix <- jamba::rbindList(
                        jamba::rmNULL(strsplit(as.character(ix), split),
                           nullValue=""));
                     colnames(ix) <- jamba::makeNames(rep(i, ncol(ix)));
                  }
               }
               ix;
            })
         )
      );
   }
   # updated to ignore intermediate and intermediate_source
   xnames <- setdiff(colnames(x),
      c(intermediate, intermediate_source));

   if (length(try_list) > 0) {
      # 30mar2021: port from "found" to intermediate
      #if (!"found" %in% colnames(x)) {
      #   x[["found"]] <- rep("", nrow(x));
      #}
      if (!intermediate %in% colnames(x)) {
         x[[intermediate]] <- rep("", nrow(x));
      }
      # 30mar2021: port from "found_source" to intermediate_source
      #if (!"found_source" %in% colnames(x)) {
      #   x[["found_source"]] <- rep("", nrow(x));
      #   if ("first_try" %in% handle_multiple) {
      #      x[["found_try"]] <- rep(TRUE, nrow(x));
      #   }
      #}
      if (!intermediate_source %in% colnames(x)) {
         x[[intermediate_source]] <- rep("", nrow(x));
         if ("first_try" %in% handle_multiple) {
            x[["found_try"]] <- rep(TRUE, nrow(x));
         }
      }
   }

   for (itry in try_list) {
      for (iann in ann_lib) {
         if (verbose) {
            jamba::printDebug("freshenGenes(): ",
               "iann:'",
               iann, "'");
         }
         if ("character" %in% class(itry)) {
            itryname <- paste0(gsub("[.]db$", "", iann), itry);
            if (verbose) {
               jamba::printDebug("freshenGenes(): ",
                  "itryname:'",
                  itryname, "'");
            }
            ienv <- get_anno_db(itryname,
               verbose=verbose,
               ...);
            itryname <- attr(ienv, "annoname");
         } else {
            if (verbose) {
               jamba::printDebug("freshenGenes(): ",
                  "Using ann_lib entry as-is");
            }
            ienv <- get_anno_db(itry,
               verbose=verbose,
               ...);
            itryname <- attr(ienv, "annoname");
         }
         if (length(ienv) == 0) {
            if (verbose) {
               jamba::printDebug("freshenGenes(): ",
                  "   Skipping", fgText=c("darkorange1", "red"));
            }
            next;
         }
         
         for (iname in xnames) {
            if (verbose) {
               jamba::printDebug("freshenGenes(): ",
                  "   iname:",
                  iname);
            }
            #ifound <- x[["found"]];
            #ifound_source <- x[["found_source"]];
            ifound <- x[[intermediate]];
            ifound_source <- x[[intermediate_source]];
            if ("first_hit" %in% handle_multiple) {
               ## input must have characters
               ## must have no found result
               #ido <- (nchar(jamba::rmNA(naValue="", ifound)) == 0 &
               #      nchar(jamba::rmNA(naValue="", x[[iname]])) > 0);
               ido <- (genejam::is_empty(ifound) &
                     !genejam::is_empty(x[[iname]]))
            } else if ("first_try" %in% handle_multiple) {
               ## input must have characters
               ## previous try must have no result
               ido <- (x[["found_try"]] &
                     !genejam::is_empty(x[[iname]]));
                     #nchar(jamba::rmNA(naValue="", x[[iname]])) > 0);
            } else {
               ## input must have characters, and not be empty
               #ido <- (nchar(jamba::rmNA(naValue="", x[[iname]])) > 0);
               ido <- !genejam::is_empty(x[[iname]]);
            }
            ix <- x[[iname]][ido];
            if (length(ix) == 0) {
               if (verbose) {
                  jamba::printDebug("freshenGenes(): ",
                     "      Skipping because ", "0", " entries to query.",
                     fgText=c("darkorange1", "red"));
               }
               next;
            }
            if (verbose) {
               jamba::printDebug("freshenGenes(): ",
                  "      Querying ",
                  jamba::formatInt(length(ix)),
                  " entries.",
                  fgText=c("darkorange1", "aquamarine3"));
            }
            ixu <- jamba::rmNA(as.character(unique(ix)));
            ivals_l <- AnnotationDbi::mget(ixu,
               ienv,
               ifnotfound=NA);

            ## Data cleaning step to protect values which have delimiters
            if (protect_inline_sep && jamba::igrepHas(sep, unlist(ivals_l))) {
               ## convert to dummy '!:!'
               if (verbose) {
                  jamba::printDebug("freshenGenes(): ",
                     "Converted intermediate '", 
                     sep, 
                     "' to '", 
                     "!:!", 
                     "'");
               }
               ivals_l <- lgsub(sep, "!:!", ivals_l);
            }
            
            ivals <- jamba::cPaste(ivals_l,
               sep=sep,
               na.rm=TRUE);
            
            names(ivals) <- ixu;
            ivals <- ivals[!is.na(ivals)];
            ixnew <- ivals[match(ix, names(ivals))];
            ixdo <- (nchar(ixnew) > 0);
            
            iname_tryname <- itryname;
            if (any(c("first_try", "all") %in% handle_multiple)) {
               # fix slight bug in handling ifound_source=""
               # that would create output like ",source_name"
               ifound_source[ido][ixdo] <- ifelse(
                  nchar(ifound[ido][ixdo]) == 0,
                  iname_tryname,
                  ifelse(nchar(ifound_source[ido][ixdo]) == 0,
                     iname_tryname,
                     paste0(ifound_source[ido][ixdo], sep, iname_tryname)));
               ifound[ido][ixdo] <- ifelse(
                  nchar(ifound[ido][ixdo]) == 0,
                  ixnew[ixdo],
                  paste0(ifound[ido][ixdo], sep, ixnew[ixdo]));
            } else {
               ifound[ido][ixdo] <- ixnew[ixdo];
               ifound_source[ido][ixdo] <- iname_tryname;
            }
            x[[intermediate]] <- ifound;
            x[[intermediate_source]] <- ifound_source;
         }
         if ("first_try" %in% handle_multiple) {
            isempty <- genejam::is_empty(x[[intermediate]]);
            x[["found_try"]][!isempty] <- FALSE;
         }
      }
   }
   
   ###################################
   ## Remove found_try column
   if ("first_try" %in% handle_multiple) {
      x <- x[,setdiff(colnames(x), "found_try"),drop=FALSE];
   }
   
   ###################################
   ## Optionally remove found_source
   if (!include_source) {
      x <- x[,setdiff(colnames(x), intermediate_source),drop=FALSE];
   }
   
   ###################################
   ## Make intermediate values unique
   ## also this step sorts delimited values
   if (intermediate %in% colnames(x) && any(nchar(x[[intermediate]]) > 0)) {
      xfoundu <- unique(x[[intermediate]]);
      if (length(split) > 0 && nchar(split) > 0) {
         xfounduv <- jamba::cPasteSU(strsplit(xfoundu, split),
            sep=sep);
      } else {
         xfounduv <- xfoundu;
      }
      x[[intermediate]] <- xfounduv[match(x[[intermediate]], xfoundu)];
   }
   
   ## Revert protected sep values
   if (protect_inline_sep && jamba::igrepHas("!:!", x[[intermediate]])) {
      ## convert from dummy '!:!' to sep
      if (verbose) {
         jamba::printDebug("freshenGenes(): ",
            "Converted intermediate '", 
            "!:!", 
            "' back to '", 
            sep, 
            "'");
      }
      x[[intermediate]] <- gsub("!:!", sep, x[[intermediate]]);
   }
   
   ## make found_source values unique, if they are retained in output
   if (include_source && intermediate_source %in% colnames(x)) {
      xfoundsu <- unique(x[[intermediate_source]]);
      if (length(split) > 0 && nchar(split) > 0) {
         xfoundsuv <- jamba::cPasteU(strsplit(xfoundsu, split),
            sep=sep);
      } else {
         xfoundsuv <- xfoundsu;
      }
      x[[intermediate_source]] <- xfoundsuv[match(x[[intermediate_source]], xfoundsu)];
   }
   if (verbose) {
      jamba::printDebug("freshenGenes(): ",
         "head(x, 10):");
      print(head(x, 10));
   }
   
   ###################################
   ## final arrangement of columns
   ## 30mar2021: changed to stop using "found"
   ## which required renaming to "intermediate"
   ## and use intermediate directly
   if (length(final) > 0) {
      #xnames <- colnames(x);
      #xnames <- jamba::makeNames(
      #   gsub("^found", 
      #      "intermediate",
      #      xnames),
      #   renameFirst=FALSE);
      #colnames(x) <- xnames;
      
      ## LOC# recovery for entries that have no intermediate
      # 30mar2021 this change in isempty should be slightly faster
      #isempty <- (nchar(jamba::rmNA(naValue="", x[[intermediate]])) == 0);
      #isempty <- (is.na(x[[intermediate]]) | nchar(x[[intermediate]]) == 0);
      isempty <- genejam::is_empty(x[[intermediate]]);
      if (any(isempty) && length(xnames) > 0) {
         xnames1 <- head(xnames, 1);
         isloc <- grepl("^LOC[0-9]+$", x[[xnames1]][isempty]);
         if (any(isloc)) {
            if (verbose) {
               jamba::printDebug("freshenGenes(): ",
                  "Converting ",
                  jamba::formatInt(sum(isloc)),
                  " entries with format ",
                  "'LOC1234567'",
                  " and no intermediate, to: ",
                  "'1234567'", " format.");
            }
            # replace "LOC211052" with "211052"
            x[[intermediate]][isempty][isloc] <- gsub("^LOC",
               "",
               x[[xnames1]][isempty][isloc]);
         }
      }
      if (verbose) {
         jamba::printDebug("freshenGenes(): ",
            "Processing final data.frame, head(x, 10):");
         print(head(x, 10));
      }
      for (i in final) {
         if (verbose) {
            jamba::printDebug("freshenGenes(): ",
               "Processing final:", i);
         }
         # note intermediate=i and final=NULL
         # means the return data will have colname i
         x1 <- freshenGenes(x[[intermediate]],
            sep=sep,
            split=sep,
            handle_multiple=handle_multiple,
            ann_lib=ann_lib,
            try_list=i,
            final=NULL,
            include_source=FALSE,
            intermediate=i,
            verbose=verbose > 1,
            ...);
         x[[i]] <- x1[[i]];
         # This step should not be necessary for "final"
         #if (include_source) {
         #   i_source <- paste0(i, "_source");
         #   x[[i_source]] <- x1[[i_source]];
         #}
      }
      if ("original" %in% empty_rule) {
         ifinal <- head(final, 1);
         # 30mar2021 this change in isempty should be slightly faster
         #isempty <- (nchar(jamba::rmNA(naValue="", x[[ifinal]])) == 0);
         #isempty <- (is.na(x[[intermediate]]) | nchar(x[[intermediate]]) == 0);
         isempty <- genejam::is_empty(x[[intermediate]]);
         x[[ifinal]][isempty] <- x[[1]][isempty];
      } else if ("na" %in% empty_rule) {
         ifinal <- head(final, 1);
         # 30mar2021 this change in isempty should be slightly faster
         #isempty <- (nchar(jamba::rmNA(naValue="", x[[ifinal]])) == 0);
         #isempty <- (is.na(x[[intermediate]]) | nchar(x[[intermediate]]) == 0);
         isempty <- genejam::is_empty(x[[intermediate]]);
         x[[ifinal]][isempty] <- NA;
      }
   }
   
   return(x);
}

#' Get annotation database or environment
#' 
#' Get annotation database or environment
#' 
#' This function is a simple wrapper function that takes either an
#' annotation data name, for example from the `AnnotationDbi` package,
#' or an annotation object, and returns the annotation object.
#' 
#' In the event the annotation object must be derived using
#' `AnnotationDbi::revmap()`, then that process is performed, and the
#' reverse mapped annotation object is returned.
#' 
#' @family genejam
#' 
#' @param x character name of an annotation object, or an annotation
#'    object itself.
#' @param revmap_suffix character string indicting the expected suffix
#'    that can be used to create reverse-mapped annotation data, for
#'    example the suffix `"2EG"` is used to indicate that annotation
#'    returns Entrez gene. When annotation does not contain this suffix,
#'    the annotation is reverse-mapped using `AnnotationDbi::revmap()`.
#' @param ignore.case logical indicating whether to return an environment
#'    after converting all keys to lowercase, which is one implementation
#'    choice to provide case-insensitive output from `mget()`. In order
#'    to fulfill the potential, the subsequent `mget()` must also
#'    use `tolower(x)` on the input character vector. Note that this
#'    option is currently fairly slow, and uses more memory while the
#'    environment is loaded.
#' @param verbose logical indicating whetheer to print verbose output.
#' @param ... additional arguments are ignored.
#' 
#' @export
get_anno_db <- function
(x,
 revmap_suffix=c("2EG", "2ENTREZID", "2NAME"),
 ignore.case=FALSE,
 verbose=FALSE,
 ...)
{
   #
   if ("character" %in% class(x)) {
      if (verbose) {
         jamba::printDebug("get_anno_db(): ",
            x);
      }
      flip_itryname <- function(x, revmap_suffix, verbose=FALSE) {
         ## This function tests if the reciprocal name exists,
         ## and if so it returns that name.
         ## Otherwise it returns NULL.
         itryname <- NULL;
         if (length(revmap_suffix) > 0 && any(nchar(revmap_suffix) > 0)) {
            revmap_suffix <- revmap_suffix[nchar(revmap_suffix) > 0];
            revmap_anygrep <- paste0("(",
               jamba::cPaste(revmap_suffix, 
                  sep="|"), 
               ")$");
            if (jamba::igrepHas(revmap_anygrep, x)) {
               ## one of the revmap extensions exists as a suffix, remove it
               itryname <- gsub(revmap_anygrep, "", x);
               if (!better_exists(itryname)) {
                  itryname <- NULL;
               }
               return(itryname);
            }
            for (revmap_suffixi in revmap_suffix) {
               itryname <- paste0(x, revmap_suffixi);
               if (verbose) {
                  jamba::printDebug("flip_itryname(): ",
                     "itryname:",
                     itryname);
               }
               if (better_exists(itryname)) {
                  return(itryname);
               }
               itryname <- NULL;
            }
         }
         return(itryname);
      }
      if (better_exists(x)) {
         ## If the name exists, return it directly
         itry <- better_get(x);
      } else {
         ## If the name does not exist, test for the reciprocal name
         reciprocal_x <- flip_itryname(x,
            revmap_suffix=revmap_suffix,
            verbose=verbose);
         if (length(reciprocal_x) > 0) {
            itry <- AnnotationDbi::revmap(better_get(reciprocal_x));
         } else {
            if (verbose) {
               jamba::printDebug("get_anno_db(): ",
                  "Not found on the search path.");
            }
            return(NULL);
         }
      }
      attr(itry, "annoname") <- x;
   } else {
      itrynames <- jamba::provigrep(c("objTarget", "objName"),
         slotNames(x));
      itryname <- paste(
         unlist(lapply(itrynames, function(i){
            slot(x, i)
         })),
         collapse=".");
      attr(x, "annoname") <- itryname;
      itry <- x;
   }
   if (ignore.case) {
      itry_l <- as.list(itry);
      names(itry_l) <- tolower(names(itry_l));
      itry <- as.environment(itry_l);
      rm(itry_l);
   }
   return(itry);
}

#' Better exists()
#' 
#' Better exists()
#' 
#' @family jam utility functions
#' 
#' This function is a lightweight enhancement of `base::exists()`
#' that accepts a package prefix in the object name, for example
#' `"base::exists"`.
#' 
#' This function recognizes a package prefix `"packagename::"`
#' and uses it to determine the correct value for argument `where`.
#' If the package is not present in `search()` using the
#' form `"package:packagename"` then an error is thrown.
#' Otherwise if the package is on the search path, this
#' function simply calls `base::exists(x, where=posnum, ...)`.
#' 
#' If the input `x` does not contain a package prefix, then
#' this function simply calls `base::exists(x)`
#' for the default behavior.
#' 
#' @return logical vector with length equal to `length(x)`. Note that
#'    when a package prefix is supplied, when the package is not on
#'    the search path this function returns `FALSE`, and does not
#'    throw an error.
#' 
#' @param x character vector length 1 or more, indicating the
#'    object names, with or without package prefix.
#' @param where,mode,inherits arguments passed to `base::exists()`.
#'    Note that arguments `envir` and `frame` use the defaults,
#'    and therefore may not be compatible with using input `x`
#'    with more than one value.
#' @param ... additional arguments are passed to `base::exists()`.
#' 
#' #' @examples
#' exists("exists", where="package:base")
#' 
#' better_exists("base::exists")
#' 
#' @export
better_exists <- function
(x,
 where=-1,
 mode="any",
 inherits=TRUE,
 verbose=FALSE,
 ...)
{
   ##
   if (length(x) > 1) {
      where <- rep(where, length.out=length(x));
      mode <- rep(mode, length.out=length(x));
      inherits <- rep(inherits, length.out=length(x));
      be_result <- sapply(seq_along(x), function(xi){
         better_exists(x[[xi]],
            where=where[[xi]],
            mode=mode[[xi]],
            inherits=inherits[[xi]],
            ...);
      });
      return(be_result);
   }
   if (jamba::igrepHas("^.+::.+$", x)) {
      xpackage <- paste0("package:", gsub("^(.+)::(.+)$", "\\1", x));
      xbase <- gsub("^(.+)::(.+)$", "\\2", x);
      xpos <- jamba::rmNA(match(xpackage, search()));
      if (verbose) {
         jamba::printDebug("better_exists(): ",
            "xpackage:", xpackage);
         jamba::printDebug("better_exists(): ",
            "xbase:", xbase);
         jamba::printDebug("better_exists(): ",
            "xpos:", xpos);
      }
      if (length(xpos) == 0) {
         return(FALSE);
         stop(paste0("better_exists(): There is no package called '", 
            xpackage, 
            "' on the search() list."));
      }
      base::exists(xbase,
         where=xpos, 
         mode=mode,
         inherits=inherits,
         ...);
   } else {
      base::exists(x,
         where=where,
         mode=mode,
         inherits=inherits,
         ...);
   }
}

#' Better get()
#' 
#' Better get()
#' 
#' @family jam utility functions
#' 
#' This function is a lightweight enhancement of `base::get()`
#' that accepts a package prefix in the object name, for example
#' `"base::exists"`.
#' 
#' This function recognizes a package prefix `"packagename::"`
#' and uses it to determine the correct value for argument `where`.
#' If the package is not present in `search()` using the
#' form `"package:packagename"` then an error is thrown.
#' Otherwise if the package is on the search path, this
#' function simply calls `base::get(x, pos=posnum, ...)`.
#' 
#' If the input `x` does not contain a package prefix, then
#' this function simply calls `base::get()`
#' for the default behavior.
#' 
#' @return the R object found. If not object is found an error results.
#' 
#' @param character string with an R object name, with or without
#'    an R package prefix. Note that only one value is recognized.
#' @param pos,mode,inherits arguments passed to `base::get()`. Note
#'    that argument `envir` is not passed to `base::get()`, since
#'    it is typically defined dynamically by that function.
#' @param ... additional arguments are passed to `base::get()`. The
#'    only additional argument is `envir` which is not recommended,
#'    but is here for compatibility with the base functionality.
#' 
#' @examples 
#' get("get", pos="package:base")
#' 
#' better_get("base::get")
#' 
#' @export
better_get <- function
(x,
 pos=-1L,
 mode="any",
 inherits=TRUE,
 ...)
{
   if (jamba::igrepHas("^.+::.+$", x)) {
      xpackage <- paste0("package:", gsub("^(.+)::(.+)$", "\\1", x));
      xbase <- gsub("^(.+)::(.+)$", "\\2", x);
      xpos <- match(xpackage, search());
      if (length(xpos) == 0) {
         stop(paste0("better_get(): There is no package called '", 
            xpackage, 
            "' on the search() list."));
      }
      base::get(xbase,
         pos=xpos, 
         mode=mode,
         inherits=inherits,
         ...);
   } else {
      base::get(x,
         pos=pos,
         mode=mode,
         inherits=inherits,
         ...);
   }
}

#' Pattern replacement in a list of character vectors
#' 
#' Pattern replacement in a list of character vectors
#' 
#' This function is a simple wrapper around `base::gsub()` except
#' it operates on a list.
#' 
#' Note that this function assumes the input data contains vectors
#' and not embedded list objects.
#' 
#' @family jam list functions
#' 
#' @param pattern,replacement,ignore.case,perl,fixed,useBytes all
#'    arguments are passed to `base::gsub()` after `x` is converted
#'    to a character vector.
#' @param x `list` object that contains character vectors.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' 
#' x <- list(a=c("A", "B"), b=c("C,D"));
#' lgsub(",", "!:!", x)
#' 
#' @export
lgsub <- function
(pattern,
 replacement,
 x,
 ignore.case=FALSE,
 perl=FALSE,
 fixed=FALSE,
 useBytes=FALSE,
 ...)
{
   ## Expand x
   xlen <- lengths(x);
   xsplit <- rep(seq_along(x), xlen);
   xexp <- unname(unlist(x));
   xexp_new <- gsub(pattern=pattern,
      replacement=replacement,
      x=xexp,
      perl=perl,
      fixed=fixed,
      useBytes=useBytes);
   x_out <- split(xexp_new, xsplit);
   if (length(names(x)) > 0) {
      names(x_out) <- names(x);
   }
   return(x_out);
}

#' Case-insensitive mget()
#' 
#' Case-insensitive mget()
#' 
#' This function is a lightweight wrapper around `base::mget()`
#' (and generics) that intends to allow case-insensitive matching.
#' It does so by converting all keys to lowercase, matching
#' lowercase input to these lowercase keys, then using the original
#' keys in native `base::mget()`.
#' 
#' One small change from `base::mget()` is the default
#' argument `ifnotfound=NA`.
#' 
#' This function secretly runs `mget()` using the unique lowercase
#' input values `x`, to reduce the number of queries. This implementation
#' is designed to help with extremely long and potentially highly duplicated
#' input values in `x`, in which case the change greatly reduces the time
#' to return results.
#' 
#' Note: This function returns the first matching lowercase
#' key, with the direct assumption that keys will not be duplicated
#' after converting to lowercase. Should this assumption become a
#' problem, please provide feedback and we will change the method
#' accordingly.
#' 
#' Note: For unknown reasons, the R method dispatch was not
#' behaving properly for objects of class `"AnnDbBimap"`, presumably
#' because the generic functions `AnnotationDbi::ls()` and
#' `AnnotationDbi::mget()` were written for class `"Bimap"`.
#' So when the input `envir` class contains `"Bimap"` the
#' direct function `AnnotationDbi::ls()` or
#' `AnnotationDbi::mget()` is called, otherwise the generic
#' `ls()` or `mget()` is called.
#' 
#' @return named `list` of objects found, or `NA` for objects
#'    that are not found.
#' 
#' @family jam utility functions
#' 
#' @param x character vector of object names.
#' @param envir,mode,ifnotfound,inherits arguments are passed to
#'    `base:mget()`.
#' 
#' @export
imget <- function
(x,
 envir=as.environment(-1L),
 mode="any",
 ifnotfound=NA,
 inherits=FALSE,
 verbose=FALSE,
 ...)
{
   ##
   if (jamba::igrepHas("Bimap", class(envir))) {
      keys <- AnnotationDbi::ls(envir);
   } else {
      keys <- ls(envir);
   }
   xl <- tolower(x);
   xlu <- unique(xl);
   xmatch <- match(xl, xlu);
   
   ## Prepare an empty output list
   valuesu <- as.list(rep(NA, length(xlu)));
   names(valuesu) <- xlu;
   
   ## Match unique lowercase input to lowercase keys
   keymatch <- match(xlu, tolower(keys));
   keysfound <- !is.na(keymatch);
   if (any(keysfound)) {
      if (jamba::igrepHas("Bimap", class(envir))) {
         valuesfound <- AnnotationDbi::mget(keys[keymatch[keysfound]],
            envir, 
            ifnotfound=ifnotfound,
            inherits=inherits);
      } else {
         valuesfound <- mget(keys[keymatch[keysfound]],
            envir, 
            ifnotfound=ifnotfound,
            inherits=inherits);
      }
      valuesu[keysfound] <- valuesfound;
   }
   values <- valuesu[xmatch];
   values;
}
