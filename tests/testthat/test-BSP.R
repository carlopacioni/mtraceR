library(mtraceR)
context("test BSP")


test_that("test skyline values", {
  # dir
  data.path <- system.file("extdata", package="mtraceR")
  
  test_skyline <- BSP (dir.in=data.path, skylinefile="skylinefile", 
                       dir.out=tempdir(), all.loci=TRUE, overall=TRUE, params=1, 
                       gen=1, mu=1, save2disk=FALSE)
  expect_equal(dim(test_skyline[[1]]), c(1443,   12))
  
  expect_equal(test_skyline[[1]][1, Parameter.value], 60.38454)
  expect_equal(test_skyline[[1]][1443, Lower], 2630.3637)
  
})