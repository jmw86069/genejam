

#' Freshen gene annotations using Bioconductor annotation data
#' 
#' Freshen gene annotations using Bioconductor annotation data
#' 
#' This function takes a vector or `data.frame` of gene symbols,
#' and uses Bioconductor annotation methods to find the most current
#' official gene symbol.
#' 
#' The function can also be customized to return
#' additional annotation, or custom annotation.
#' 
#' The annotation process is intended to run in two steps:
#' 
#' 1. Convert the input gene symbol to the recognized Entrez gene ID,
#' usually using something like `"org.Hs.egSYMBOL2EG"` which returns
#' a character string with Entrez gene identifier. When the
#' gene symbol is not recognized, the backup plan is to use gene
#' aliases, accession numbers, or other relevant annotations to
#' find the current official gene. However, these steps are not
#' as trusted as using the `"SYMBOL2EG"` approach, so that step
#' is attempted first.
#' 2. Convert the Entrez gene identifier to the official Entrez gene
#' symbol, usually using something like `"org.Hs.egSYMBOL"`. This
#' step is also a convenient time to include things like the
#' descriptive gene name, which is provided by `"org.Hs.egGENENAME"`.
#' 
#' @return `data.frame` with one or more columns indicating the input
#' data, then a column `"intermediate"` containing the Entrez gene ID
#' that was matched, then one column for each item in `finalList`,
#' by default `"SYMBOL"`.
#' 
#' @family genejam
#' 
#' @param x character vector or `data.frame` with one or most columns
#'    containing gene symbols.
#' @param annLib character value indicating the name of the Bioconductor
#'    annotation library to use when looking up gene nomenclature.
#' @param tryList character vector indicating one or more names of
#'    annotations to use for the input gene symbols in `x`. The
#'    annotation should typically return the Entrez gene ID, usually
#'    given by `'2EG'` at the end of the name. For example `SYMBOL2EG`
#'    will be used with annLib `"org.Hs.eg.db"` to produce annotation
#'    name `"org.Hs.egSYMBOL2EG"`. Note that when the `'2EG'` form of
#'    annotation does not exist, it will be derived using
#'    `AnnotationDbi::revmap()`. For example if `"org.Hs.egALIAS"`
#'    exists, but not `"org.Hs.egALIAS2EG"`, then this function will
#'    create a reverse-mapped `"org.Hs.egALIAS2EG"` derived from
#'    `"org.Hs.egALIAS"`.
#' @param finalList character vector to use for the final conversion
#'    step. When finalList is `NULL` no conversion is performed.
#'    When `finalList` contains multiple values, each value is returned
#'    in the output. For example, `finalList=c("SYMBOL","GENENAME")` will
#'    return a column `"SYMBOL"` and a column `"GENENAME"`.
#' @param sep character value used to separate delimited values in `x`.
#'    For example when `sep=","` then comma-delimited values will be split
#'    into separate columns, and each column can be used in the gene
#'    annotation update. See `handle_multiple`.
#' @param handle_multiple character value indicating how to handle multiple
#'    values: `"first_hit"` will query each column of `x` until it finds the
#'    first possible returning match, and will ignore all subsequent possible
#'    matches for that row in `x`. For example, if one row in `x` contains
#'    multiple values, only the first match will be used. `"first_try"`
#'    will return the first match from `tryList` for all columns in `x`
#'    that contain a match. For example, if one row in `x` contains two
#'    values, the first match from `tryList` using one or both columns in
#'    `x` will be maintained. Subsequent entries in `tryList` will not be
#'    attempted for rows that already have a match. `"all"` will return all
#'    possible matches for all entries in `x` using all items in `tryList`.
#' @param empty_rule character value indicating how to handle entries which
#'    did not have a match, and are therefore empty: `"original"` will use
#'    the original entry as the output field; `"empty"` will leave the
#'    entry blank.
#' @param include_source logical indicating whether to include a column
#'    that shows the colname and source matched. For example, if column
#'    `"original_gene"` matched `"SYMBOL2EG"` in `"org.Hs.eg.db"` there
#'    will be a column `"found_source"` with value
#'    `"original_gene.org.Hs.egSYMBOL2EG"`.
#' @param verbose logical indicating whether to print verbose output.
#' 
#' @examples
#' if (suppressPackageStartupMessages(require(org.Hs.eg.db))) {
#'    freshenGenes(c("APOE", "CCN2", "CTGF"));
#'    
#'    ## Optionally show the annotation source matched
#'    freshenGenes(c("APOE", "CCN2", "CTGF"), include_source=TRUE)
#'    
#'    ## Show comma-delimited genes
#'    freshenGenes(c("APOE", "CCN2", "CTGF", "CCN2,CTGF"));
#'    
#'    ## Optionally include more than SYMBOL in the output
#'    freshenGenes(c("APOE", "CCN2", "CTGF"),
#'       finalList=c("SYMBOL", "ALIAS", "GENENAME"))
#' }
#' 
#' @export
freshenGenes <- function
(x,
 annLib=c("org.Hs.eg.db"),
 tryList=c("SYMBOL2EG", "ACCNUM2EG", "ALIAS2EG"),
 finalList=c("SYMBOL"),
 split="[, /]+",
 sep=",",
 handle_multiple=c("first_try", "first_hit", "all"),
 empty_rule=c("original", "empty"),
 include_source=FALSE,
 verbose=FALSE,
 ...)
{
   ###
   handle_multiple <- match.arg(handle_multiple);
   empty_rule <- match.arg(empty_rule);
   if (!suppressPackageStartupMessages(require(annLib, character.only=TRUE))) {
      stop("Not all packages in annLib are available.");
   }
   if (is.atomic(x)) {
      x <- data.frame(input=as.character(x),
         stringsAsFactors=FALSE);
   }
   if (length(colnames(x)) == 0) {
      colnames(x) <- jamba::makeNames(rep("input", ncol(x)));
   }
   ## Expand columns containing delimited values if necessary
   if (length(split) > 0) {
      x <- data.frame(stringsAsFactors=FALSE,
         do.call(cbind,
         lapply(jamba::nameVector(colnames(x)), function(i){
            ix <- as.character(x[[i]]);
            if (jamba::igrepHas(split, ix)) {
               ix <- jamba::rbindList(
                  jamba::rmNULL(strsplit(as.character(ix), split),
                     nullValue=""));
               colnames(ix) <- jamba::makeNames(rep(i, ncol(ix)));
            }
            ix;
         })));
   }
   xnames <- colnames(x);
   x[["found"]] <- rep("", nrow(x));
   x[["found_source"]] <- rep("", nrow(x));
   if ("first_try" %in% handle_multiple) {
      x[["found_try"]] <- rep(TRUE, nrow(x));
   }
   
   ## Iterate each column to find a match
   for (itry in tryList) {
      if (verbose) {
         jamba::printDebug("itry:", itry);
      }
      if ("character" %in% class(itry)) {
         itryname <- paste0(gsub("[.]db$", "", annLib), itry);
         if (verbose) {
            jamba::printDebug(itryname);
         }
         itry <- get_anno_db(itryname);
         itryname <- attr(itry, "annoname");
      } else {
         itry <- get_anno_db(itry);
         itryname <- attr(itry, "annoname");
      }
      for (iname in xnames) {
         if (verbose) {
            jamba::printDebug("   iname:", iname);
         }
         ifound <- x[["found"]];
         ifound_source <- x[["found_source"]];
         if ("first_hit" %in% handle_multiple) {
            ## input must have characters
            ## must have no found result
            ido <- (nchar(ifound) == 0 & nchar(x[[iname]]) > 0);
         } else if ("first_try" %in% handle_multiple) {
            ## input must have characters
            ## previous try must have no result
            ido <- (x[["found_try"]] & nchar(x[[iname]]) > 0);
         } else {
            ## input must have characters
            ido <- (nchar(x[[iname]]) > 0);
         }
         ix <- x[[iname]][ido];
         ixu <- jamba::rmNA(as.character(unique(ix)));
         if (verbose) {
            jamba::printDebug("ix:", head(ix, 10));
            jamba::printDebug("ixu:", head(ixu, 10));
            jamba::printDebug("class(itry):", class(itry));
         }
         ivals_l <- AnnotationDbi::mget(ixu,
            itry,
            ifnotfound=NA);
         ivals <- jamba::cPaste(ivals_l,
            na.rm=TRUE);
         #printDebug("      ivals:", ivals);
         names(ivals) <- ixu;
         ivals <- ivals[!is.na(ivals)];
         ixnew <- ivals[match(ix, names(ivals))];
         ixdo <- (nchar(ixnew) > 0);
         #iname_tryname <- paste0(iname, ":", itryname);
         iname_tryname <- itryname;
         if (any(c("first_try","all") %in% handle_multiple)) {
            ifound_source[ido][ixdo] <- ifelse(
               nchar(ifound[ido][ixdo]) == 0,
               iname_tryname,
               paste0(ifound_source[ido][ixdo], sep, iname_tryname));
            ifound[ido][ixdo] <- ifelse(
               nchar(ifound[ido][ixdo]) == 0,
               ixnew[ixdo],
               paste0(ifound[ido][ixdo], sep, ixnew[ixdo]));
         } else {
            ifound[ido][ixdo] <- ixnew[ixdo];
            ifound_source[ido][ixdo] <- iname_tryname;
         }
         x[["found"]] <- ifound;
         x[["found_source"]] <- ifound_source;
      }
      if ("first_try" %in% handle_multiple) {
         isnonempty <- (nchar(jamba::rmNA(naValue="", x[["found"]])) > 0);
         x[["found_try"]][isnonempty] <- FALSE;
      }
   }
   ## Remove found_try column
   if ("first_try" %in% handle_multiple) {
      x <- x[,setdiff(colnames(x), "found_try"),drop=FALSE];
   }
   ## Remove found_source
   if (!include_source) {
      x <- x[,setdiff(colnames(x), "found_source"),drop=FALSE];
   }
   
   ## Make values unique
   xfoundu <- unique(x[["found"]]);
   xfounduv <- jamba::cPasteSU(strsplit(xfoundu, sep));
   x[["found"]] <- xfounduv[match(x[["found"]], xfoundu)];
   if (include_source) {
      xfoundsu <- unique(x[["found_source"]]);
      xfoundsuv <- jamba::cPasteU(strsplit(xfoundsu, sep));
      x[["found_source"]] <- xfoundsuv[match(x[["found_source"]], xfoundsu)];
   }
   
   ###############################
   ## finalList
   if (length(finalList) > 0) {
      xnames <- colnames(x);
      xnames <- jamba::makeNames(gsub("^found", 
         "intermediate", 
         xnames));
      colnames(x) <- xnames;
      ## LOC# recovery
      isempty <- (nchar(jamba::rmNA(naValue="", x[["intermediate"]])) == 0);
      if (any(isempty)) {
         isloc <- grepl("^LOC[0-9]+$", x[[1]][isempty]);
         if (any(isloc)) {
            x[["intermediate"]][isempty][isloc] <- gsub("^LOC", "", x[[1]][isempty][isloc]);
         }
      }
      for (i in finalList) {
         if (verbose) {
            jamba::printDebug("final i:", i);
         }
         x1 <- freshenGenes(x[["intermediate"]],
            sep=sep,
            handle_multiple=handle_multiple,
            annLib=annLib,
            tryList=i,
            finalList=NULL,
            ...);
         x[[i]] <- x1[["found"]];
         if (include_source) {
            x[[paste0(i, "_source")]] <- x1[["found_source"]];
         }
      }
      if ("original" %in% empty_rule) {
         ifinal <- head(finalList, 1);
         isempty <- (nchar(jamba::rmNA(naValue="", x[[ifinal]])) == 0);
         x[[ifinal]][isempty] <- x[[1]][isempty];
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
#' @param verbose logical indicating whetheer to print verbose output.
#' @param ... additional arguments are ignored.
#' 
#' @export
get_anno_db <- function
(x,
 revmap_suffix="2EG",
 verbose=FALSE,
 ...)
{
   #
   if ("character" %in% class(x)) {
      if (verbose) {
         jamba::printDebug("get_anno_db(): ",
            x);
      }
      if (exists(x)) {
         itry <- get(x);
      } else {
         if (length(revmap_suffix) > 0 && nchar(revmap_suffix) > 0) {
            revmap_grep <- paste0(revmap_suffix, "$");
            if (jamba::igrepHas(revmap_grep, x)) {
               itryname <- gsub(revmap_grep, "", x);
            } else {
               itryname <- paste0(x, revmap_suffix);
            }
            if (verbose) {
               jamba::printDebug("get_anno_db(): ",
                  "Applying revmap_suffix:",
                  revmap_suffix,
                  ", itryname:",
                  itryname);
            }
            if (exists(itryname)) {
               itry <- AnnotationDbi::revmap(get(itryname));
            } else {
               stop("The annotation data was not found in the search path.");
            }
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
   return(itry);
}
