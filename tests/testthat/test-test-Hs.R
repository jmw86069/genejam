testthat::test_that("test org.Hs.eg.db", {
   # skip if org.Hs.eg.db is not available
   testthat::skip_if_not_installed("org.Hs.eg.db")
   
   tvalues <- c("CCL8", "ADCK1", "CCL8,ADCK1", "MCP-2", "MCP-2,ADCK1")
   ann_lib <- "org.Hs.eg.db"
   
   testthat::expect_equal(
      freshenGenes(tvalues, ann_lib=ann_lib)$SYMBOL,
      c("CCL8", "ADCK1", "ADCK1,CCL8", "CCL8", "ADCK1"))

   testthat::expect_equal(
      freshenGenes(tvalues,
         ann_lib=ann_lib,
         handle_multiple="first_hit")$SYMBOL,
      c("CCL8", "ADCK1", "CCL8", "CCL8", "ADCK1"))

   # "all" sorts values
   testthat::expect_equal(
      freshenGenes(tvalues,
         ann_lib=ann_lib,
         handle_multiple="al")$SYMBOL,
      c("CCL8", "ADCK1", "ADCK1,CCL8", "CCL8", "ADCK1,CCL8"))
   
   # Note best_each does not sort, keeps them in original order
   tdf <- freshenGenes3(tvalues,
      ann_lib=ann_lib,
      intermediate="ENTREZID",
      include_source=TRUE,
      handle_multiple="best_each");
   testthat::expect_equal(
      tdf$SYMBOL,
      c("CCL8", "ADCK1", "CCL8,ADCK1", "CCL8", "CCL8,ADCK1"))
   testthat::expect_equal(
      colnames(tdf),
      c("input", "ENTREZID", "ENTREZID_source", "SYMBOL", "GENENAME", "ALIAS"))
   
   tdf2 <- freshenGenes(tvalues,
      ann_lib=ann_lib,
      revert_split=FALSE,
      handle_multiple="first_try");

   idf <- data.frame(
      gene=c("CCL8,ADCK1", "", "", "CCL8"),
      ENTREZID=c("6355", "57143", "6355,57143", ""))
   idf2 <- freshenGenes(idf,
      include_source=TRUE,
      intermediate="ENTREZID",
      ann_lib=ann_lib,
      handle_multiple="best_each")
   testthat::expect_equal(
      idf2$SYMBOL,
      c("CCL8", "ADCK1", "ADCK1,CCL8", "CCL8"))
   testthat::expect_equal(
      idf2$ENTREZID_source,
      c("", "", "", "org.Hs.egSYMBOL2EG"))
   
})

testthat::test_that("test org.Mm.eg.db", {
   # skip if org.Hs.eg.db is not available
   testthat::skip_if_not_installed("org.Mm.eg.db")
   
   ann_lib_m <- "org.Mm.eg.db";
   testthat::expect_equal(
      freshenGenes(tvalues, ann_lib=ann_lib_m)$SYMBOL,
      c("", "", "", "Ccl8", "Ccl8"))

   testthat::expect_equal(
      freshenGenes(tvalues,
         ignore.case=TRUE,
         ann_lib=ann_lib_m)$SYMBOL,
      c("Ccl8", "Adck1", "Adck1,Ccl8", "Ccl8,Mcpt2", "Adck1"))
   
})
