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

test_that("mergeMS noise and maxpeaks arguments work",{
  
  MS2 <- matrix(c(100, 100.5, 101, 101.5,
                  10, 100, 10000, 1000), ncol = 2,
                dimnames = list(NULL, c("mz","intensity")))
  
  #maxpeaks tests
  expect_equal(mergeMS(MS2, maxpeaks = 9, noiselevel = 0),
               MS2)
  
  expect_equal(mergeMS(MS2, maxpeaks = 2, noiselevel = 0),
               MS2[3:4,])
  
  expect_equal(mergeMS(MS2, maxpeaks = 0, noiselevel = 0),
               MS2[FALSE,])
  
  #noiselevel tests
  expect_equal(mergeMS(MS2, maxpeaks = NULL, noiselevel = 0.1),
               MS2[3:4,])
  
  expect_equal(mergeMS(MS2, maxpeaks = NULL, noiselevel = 10),
               MS2[FALSE,])
})

test_that("mergeMS removeZeros argument works as expected",{
  
  MS2 <- matrix(c(100, 100.5, 101, 101.5,
                  10, 100, 10000, 1000), ncol = 2,
                dimnames = list(NULL, c("mz","intensity")))
  
  MS3 <- matrix(c(100, 100.25, 100.5,  100.7, 100.75,  101, 101.05, 
                   10,      0,   100,      0,      0,  10000,  0),
                ncol = 2,
                dimnames = list(NULL, c("mz","intensity")))
  
  #maxpeaks tests
  expect_equal(mergeMS(MS2, maxpeaks = 9, noiselevel = 0,
                       removeZeros = F),
               MS2)
  
  expect_equal(mergeMS(MS3, maxpeaks = NULL, noiselevel = 0,
                       removeZeros = TRUE),
               MS2[-nrow(MS2),])
  
  expect_equal(mergeMS(MS3, maxpeaks = NULL, noiselevel = 0,
                       removeZeros = FALSE),
               MS3)
  
  expect_equal(mergeMS(MS3, maxpeaks = NULL, noiselevel = 0, mzdiff = 0.1,
                       removeZeros = FALSE),
               MS3[c(-5,-nrow(MS3)),])
  
  
})

MS4 <- matrix(c(100, 100.25, 100.5,  100.7, 100.75,  101, 101.05, 
                10,      10,    10,     10,      10,  10,     10),
              ncol = 2,
              dimnames = list(NULL, c("mz","intensity")))

test_that("mergeMS toleranceFactor argument works as expected",{
  
  #tolerance factor tests
  expect_false(identical(mergeMS(MS4, noiselevel = 0,
                       removeZeros = F, mzdiff = 0.1,
                       toleranceFactor = 1),
                       mergeMS(MS4, noiselevel = 0,
                               removeZeros = F, mzdiff = 0.1,
                               toleranceFactor = 5)))
  
  
  expect_equal(mergeMS(MS4, noiselevel = 0,
                       removeZeros = F, mzdiff = 0.1,
                       toleranceFactor = c(0.4,5)),
               mergeMS(mergeMS(MS4, noiselevel = 0,
                               removeZeros = F, mzdiff = 0.1,
                               toleranceFactor = c(0.4)), noiselevel = 0,
                       removeZeros = F, mzdiff = 0.1,
                       toleranceFactor = c(5)))
  
})

test_that("mergeMS iterative argument works as expected",{
  #iterative
  MS5 <- matrix(c(100, 100.5, 101.5, 102,
                  1000, 2000, 1000, 2000), ncol = 2,
                dimnames = list(NULL, c("mz","intensity")))
  
  expect_equal(mergeMS(MS5, maxpeaks = NULL, mzdiff = 0.6,
                       noiselevel = 0, iterative = TRUE,
                       removeZeros = TRUE),
               mergeMS(MS5, maxpeaks = NULL, mzdiff = 0.6,
                       noiselevel = 0, iterative = FALSE,
                       removeZeros = TRUE))
               
  expect_equal(mergeMS(MS4, maxpeaks = NULL, mzdiff = 0.3,
          noiselevel = 0, iterative = FALSE,
          removeZeros = TRUE),
  matrix(c(mean(MS4[,1]), 70), ncol = 2,
                      dimnames = list(NULL, c("mz","intensity"))))
  
  expect_equal(mergeMS(MS4, maxpeaks = NULL, mzdiff = 0.3,
          noiselevel = 0, iterative = TRUE,
          removeZeros = TRUE),
  matrix(c(mean(MS4[1:2,1]), mean(MS4[3:5,1]), mean(MS4[6:7,1]),
                         20,               30,              20),
         ncol = 2, dimnames = list(NULL, c("mz","intensity"))))
  
  MS6 <- matrix(c(100, 100.25, 100.5,  102.7, 102.75,  103, 103.05, 
                  10,      10,    10,     10,      10,  10,     10),
                ncol = 2,
                dimnames = list(NULL, c("mz","intensity")))
  
  expect_equal(mergeMS(MS6, maxpeaks = NULL, mzdiff = 0.3,
                       noiselevel = 0, iterative = FALSE,
                       removeZeros = TRUE),
               matrix(c(mean(MS6[1:3,1]), mean(MS6[4:7,1]),
                                     30,               40), ncol = 2,
                      dimnames = list(NULL, c("mz","intensity"))))
  
  expect_equal(mergeMS(MS6, maxpeaks = NULL, mzdiff = 0.3,
                       noiselevel = 0, iterative = TRUE,
                       removeZeros = TRUE),
               matrix(c(mean(MS6[1:2,1]), MS6[3,1], mean(MS6[4:7,1]),
                                      20,       10,              40), ncol = 2,
                      dimnames = list(NULL, c("mz","intensity"))))
  
  
  MS7 <- matrix(c(100, 100.25, 100.5,  102.7, 102.75,  103, 103.05, 110,
                  10,      10,    10,     10,      10,  10,     10,  10),
                ncol = 2,
                dimnames = list(NULL, c("mz","intensity")))
  
  expect_equal(mergeMS(MS7, maxpeaks = NULL, mzdiff = 0.3,
                       noiselevel = 0, iterative = FALSE,
                       removeZeros = TRUE),
               matrix(c(mean(MS7[1:3,1]), mean(MS7[4:7,1]), 110,
                                      30,             40,   10), ncol = 2,
                      dimnames = list(NULL, c("mz","intensity"))))
  
  expect_equal(mergeMS(MS7, maxpeaks = NULL, mzdiff = 0.3,
                       noiselevel = 0, iterative = TRUE,
                       removeZeros = TRUE),
               matrix(c(mean(MS7[1:2,1]), MS7[3,1], mean(MS7[4:7,1]), 110,
                                      20,       10,              40,  10), 
                      ncol = 2, dimnames = list(NULL, c("mz","intensity"))))
})


test_that("findPatterns works",{
  
  MS2 <- matrix(c(100, 100.5, 101, 101.5,
                  10, 100, 10000, 1000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
  
  #maxpeaks tests
  expect_equal(findPatterns(list(MS2), patterns = list(pattern1 = c(100,100.5)), ppm = 5, mzdiff = 0),
               list(c(pattern1 = TRUE)))
  
  expect_equal(findPatterns(list(MS2), patterns = list(pattern1 = c(100,100.5),
                                                       pattern2 = c(100,101.5)),
                            ppm = 5, mzdiff = 0),
               list(c(pattern1 = TRUE,
                      pattern2 = TRUE)))
  
  expect_equal(findPatterns(list(MS2, NULL), patterns = list(pattern1 = c(100,100.5),
                                                       pattern2 = c(100,101.5)),
                            ppm = 5, mzdiff = 0),
               list(c(pattern1 = TRUE,
                      pattern2 = TRUE),
                    c(pattern1 = FALSE,
                      pattern2 = FALSE)))
  
  expect_equal(findPatterns(list(MS2, NULL), patterns = list(pattern1 = c(100,100.5),
                                                             pattern2 = c(100,103.5)),
                            ppm = 5, mzdiff = 0),
               list(c(pattern1 = TRUE,
                      pattern2 = FALSE),
                    c(pattern1 = FALSE,
                      pattern2 = FALSE)))
  })