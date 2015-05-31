library(mtraceR)
context("test mtrace")


test_that("test seq custom matrix", {
  # dir
  data.path <- system.file("extdata", package="mtraceR")
  
  # load data
  data(mtDNA_4popsSt)
  suppressMessages(test_mtDNA <- mtrace(burn.in = 0.2, save2disk = FALSE, 
                                 dir.in = data.path, 
                                 bayesallfile = "bayesallfile_mtDNA_4popsSt.gz"))
  
  data(ms_3popsSt)
  suppressMessages(test_microsats <- mtrace(burn.in = 0.4, dir.in = data.path, 
                                     save2disk = FALSE,
                                     bayesallfile = "bayesallfile_ms_3popsSt.gz"))
  
  expect_equal(test_mtDNA , mtDNA_4popsSt)
  expect_equal(test_microsats, ms_3popsSt)
})