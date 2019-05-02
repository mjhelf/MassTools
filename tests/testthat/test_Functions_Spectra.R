context("merging mass spectra")

test_that("mergeMF merges simple spectra as expected",{
  
  MS1 <- matrix(c(100, 115, 150,
                  1000, 2000, 1000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
  
  MS1merge1 <- matrix(c(110, 150,
                        3000, 1000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
 
 #merging with itself returns same spectrum if all peaks different by more than tolerance
 expect_that(MS1, equals(mergeMS(MS1)))
 expect_that(MS1, equals(mergeMS(list(MS1, MS1))))
 
 #increase tolerance to force merging
 expect_that(MS1merge1, equals(mergeMS(list(MS1), mzdiff = 16)))
 
 #same result when loading the same spectrum twice:
 expect_that(mergeMS(list(MS1), mzdiff = 16), equals(mergeMS(list(MS1, MS1), mzdiff = 16)))
 
 
  
2})