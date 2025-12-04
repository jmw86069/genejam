testthat::test_that("test org.Hs.eg.db empty_rules", {
   # skip if org.Hs.eg.db is not available
   testthat::skip_if_not_installed("org.Hs.eg.db")
   
   idf <- data.frame(
      Gene=c("APOE", "APOE", "", "TRNAV-CAC", "LOC141374963", "LOC348"),
      ENTREZID=c("", NA, "348", "", "", ""))
   outdf <- freshenGenes3(idf, include_source=TRUE, empty_rule="original")
   
   testthat::expect_equal(
      outdf$SYMBOL,
      c("APOE", "APOE", "APOE", "TRNAV-CAC", "LOC141374963", "APOE"))

   testthat::expect_equal(
      outdf$ENTREZID_source,
      c("org.Hs.egSYMBOL2EG", "", "",
         "org.Hs.egSYMBOL2EG", "", ""))

   testthat::expect_equal(
      outdf$ENTREZID,
      c("348", "348", "348",
         "107985614,107985615", "141374963", "348"))
   
})