
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
      verbose=verbose,
      ...);
}
