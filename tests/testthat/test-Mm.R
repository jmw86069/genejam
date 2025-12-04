# tests Mouse, also case-insensitive match

testthat::test_that("test org.Mm.eg.db", {
   tvalues <- c("CCL8", "ADCK1", "CCL8,ADCK1", "MCP-2", "MCP-2,ADCK1")
   
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
