context("merging mass spectra")

test_that("mergeMS merges simple spectra as expected",{
  
  MS1 <- matrix(c(100, 115, 150,
                  1000, 2000, 1000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
  
  MS1merge1 <- matrix(c(110, 150,
                        3000, 1000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
  
  MS1merge2 <- matrix(c(120,
                        4000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
 
 #merging with itself returns same spectrum if all peaks different by more than tolerance
 expect_that(MS1, equals(mergeMS(MS1)))
 expect_that(MS1, equals(mergeMS(list(MS1, MS1))))
 
 #increase tolerance to force merging
 expect_that(MS1merge1, equals(mergeMS(list(MS1), mzdiff = 16)))
 
 #same result when loading the same spectrum twice:
 expect_that(mergeMS(list(MS1), mzdiff = 16), equals(mergeMS(list(MS1, MS1), mzdiff = 16)))
 
 #same result when combining into a single peak iteratively:
 expect_that(MS1merge2, equals(mergeMS(list(MS1), mzdiff = 50)))
 expect_that(MS1merge2, equals(mergeMS(list(MS1merge1), mzdiff = 50)))
  
})

test_that("mergeMS splits peaks as expected",{
  
  MS2 <- matrix(c(100, 100.5, 101, 101.5,
                  1000, 1000, 1000, 1000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
  
  MS2merge1 <- matrix(c(100.25, 101.25,
                        2000, 2000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
  
  #merging with splitting effect from iterative peak processing
  expect_that(MS2merge1, equals(mergeMS(MS2, mzdiff = 0.5)))

  #ppm selection works, too:
  expect_that(MS2, equals(mergeMS(MS2, mzdiff = 0, ppm = 0.4*1e4)))
  expect_that(MS2merge1, equals(mergeMS(MS2, mzdiff = 0, ppm = 0.5*1e4)))
  
  })

test_that("mergeMS counts combined peaks",{
  
  MS1 <- matrix(c(100, 115, 150,
                  1000, 2000, 1000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
  
  MS1merge1 <- matrix(c(110, 150,
                        3000, 1000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
  
  MS1merge2 <- matrix(c(120,
                        4000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
  
  #merging with itself returns same spectrum if all peaks different by more than tolerance
  expect_that(c(1,1,1), equals(mergeMS(MS1, count = T)[,3]))
  expect_that(c(2,2,2), equals(mergeMS(list(MS1, MS1), count = T)[,3]))
  
  #increase tolerance to force merging
  expect_that(c(2,1), equals(mergeMS(list(MS1), mzdiff = 16, count = T)[,3]))
  
  #double the result when loading the same spectrum twice:
  expect_that(c(4,2), equals(mergeMS(list(MS1, MS1), mzdiff = 16, count = T)[,3]))
  
  #result when combining into a single peak iteratively:
  expect_that(c(count = 3), equals(mergeMS(list(MS1), mzdiff = 50, count = T)[,3]))
  expect_that(c(count = 2), equals(mergeMS(list(MS1merge1), mzdiff = 50, count = T)[,3]))
  
})
