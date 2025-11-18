
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
#' direct function `AnnotationDbi::keys()` is called, and if it
#' fails for some reason, `AnnotationDbi::ls()` is called,
#' thes `AnnotationDbi::mget()` is called, otherwise the generic
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
      keys <- tryCatch({
         AnnotationDbi::keys(envir,
            ...);
      }, error=function(e){
         AnnotationDbi::ls(envir);
      });
   } else {
      keys <- ls(envir);
   }
   xl <- tolower(x);
   xlu <- unique(xl);
   xmatch <- match(xl, xlu);
   
   ## Prepare an empty output list, not using NA
   valuesu <- rep(list(character(0)), length(xlu));
   names(valuesu) <- xlu;
   
   ## Subset relevant keys
   subkeys <- keys[tolower(keys) %in% xlu];
   
   ## Check for duplicated case-insensitive keys
   if (any(duplicated(tolower(subkeys)))) {
      if (verbose) {
         jamba::printDebug("imget(): ",
            "Detected ",
            jamba::formatInt(sum(duplicated(tolower(subkeys)))),
            " duplicate case-insensitive keys, mitigating.");
      }
      lcdupekeys <- unique(tolower(subkeys)[duplicated(tolower(subkeys))])
      alldupekeys <- subkeys[tolower(subkeys) %in% lcdupekeys]
      alldupekeys_list <- split(alldupekeys, tolower(alldupekeys))
      #
      ivals <- AnnotationDbi::mget(alldupekeys, envir)
      dupekey_list <- split(unlist(unname(ivals)),
         rep(tolower(names(ivals)), lengths(ivals)))
      # sorts
      dupekey_list <- jamba::mixedSorts(dupekey_list)
      imatch <- match(names(dupekey_list), xlu)
      valuesu[imatch] <- dupekey_list;
   }

   ## Match unique lowercase input to lowercase keys
   keymatch <- match(xlu, tolower(keys));

   # if values were already assigned, set them NA here to avoid re-assigning
   keymatch[lengths(valuesu) > 0] <- NA;
   
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
