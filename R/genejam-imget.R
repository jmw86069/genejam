
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
