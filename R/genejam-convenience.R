
#' Freshen gene annotations using Bioconductor annotation data
#' 
#' Freshen gene annotations using Bioconductor annotation data
#' 
#' This function is a convenient extension of `freshenGenes()`
#' that adds `GENENAME` to the default value for
#' `final=c("SYMBOL", "GENENAME")`.
#' It therefore returns two (`2`) annotation columns by default,
#' the gene symbol, and the long gene name.
#' 
#' @family genejam
#' 
#' @inheritParams freshenGenes
#' 
#' @export
freshenGenes2 <- function
(x,
 ann_lib=c("","org.Hs.eg.db"),
 try_list=c("SYMBOL2EG", "ACCNUM2EG", "ALIAS2EG"),
 final=c("SYMBOL", "GENENAME"),
 split="[ ]*[,/;]+[ ]*",
 sep=",",
 handle_multiple=c("first_try", "first_hit", "all", "best_each"),
 empty_rule=c("empty", "original", "na"),
 include_source=FALSE,
 protect_inline_sep=TRUE,
 intermediate="ENTREZID",
 verbose=FALSE,
 ...)
{
   freshenGenes(x,
      ann_lib=ann_lib,
      try_list=try_list,
      final=final,
      split=split,
      sep=sep,
      handle_multiple=handle_multiple,
      empty_rule=empty_rule,
      include_source=include_source,
      protect_inline_sep=protect_inline_sep,
      intermediate=intermediate,
      verbose=verbose,
      ...);
}


#' Freshen gene annotations using Bioconductor annotation data
#' 
#' Freshen gene annotations using Bioconductor annotation data
#' 
#' This function is a convenient extension of `freshenGenes()`
#' that adds `GENENAME` and `ALIAS` to the default value for
#' `final=c("SYMBOL", "GENENAME", "ALIAS")`.
#' It therefore returns three (`3`) annotation columns by default,
#' the gene symbol, the long gene name, and the common gene aliases.
#' The gene aliases often includes numerous previous gene symbols
#' attributed to the gene.
#' 
#' @family genejam
#' 
#' @inheritParams freshenGenes
#' 
#' @export
freshenGenes3 <- function
(x,
 ann_lib=c("","org.Hs.eg.db"),
 try_list=c("SYMBOL2EG", "ACCNUM2EG", "ALIAS2EG"),
 final=c("SYMBOL", "GENENAME", "ALIAS"),
 split="[ ]*[,/;]+[ ]*",
 sep=",",
 handle_multiple=c("first_try", "first_hit", "all", "best_each"),
 empty_rule=c("empty", "original", "na"),
 include_source=FALSE,
 protect_inline_sep=TRUE,
 intermediate="ENTREZID",
 verbose=FALSE,
 ...)
{
   freshenGenes(x,
      ann_lib=ann_lib,
      try_list=try_list,
      final=final,
      split=split,
      sep=sep,
      handle_multiple=handle_multiple,
      empty_rule=empty_rule,
      include_source=include_source,
      protect_inline_sep=protect_inline_sep,
      intermediate=intermediate,
      verbose=verbose,
      ...);
}


#' Test if vector elements are empty
#' 
#' Test if vector elements are empty
#' 
#' This function simply checks if values in a vector are
#' `NA` or `""` with `nchar(x) == 0`.
#' 
#' For `factor` input, the values are coerced with `as.character()`.
#' It might be slightly faster to test factor levels then to
#' apply to the full vector.
#' 
#' Todo: Make this function work with `list` input, so it
#' requires all elements to be `is_empty()`.
#' 
#' @family genejam
#' 
#' @param x `vector` that may contain `NA` values.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' x1 <- c("A", "", NA, "B,C");
#' is_empty(x1)
#' 
#' is_empty(factor(x1))
#' 
#' @export
is_empty <- function
(x,
 ...)
{
   if (!is.atomic(x)) {
      stop("Input must be atomic vector.")
   }
   if (is.factor(x)) {
      (is.na(x) | nchar(as.character(x)) == 0)
   } else if (is.numeric(x)) {
      is.na(x)
   } else {
      (is.na(x) | nchar(x) == 0)
   }
}
