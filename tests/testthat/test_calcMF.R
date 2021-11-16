context("Calculating and filterng molecula formulas")

test_that("calcMF correctly passes through Rdisop::decomposeMass results",{
  
  mtres <- calcMF(mz = 200.000659,
                  z = 1,
                  ppm = 5,
                  top = NULL,
                  elements = Rdisop::initializeCHNOPS(),
                  maxCounts = TRUE,
                  SeniorRule = TRUE,
                  HCratio = TRUE,
                  moreRatios = TRUE,
                  elementHeuristic = TRUE,
                  Filters = list(),
                  summarize = F,
                  BPPARAM = NULL)
  
  rdires <- Rdisop::decomposeMass(200.000659+ 1*5.48579909070e-4,
                                  z = 1,
                                  ppm = 5,
                                  mzabs = 0,
                                  elements = Rdisop::initializeCHNOPS())
  
  rdorder <- order(abs(rdires$exactmass - 200.000659 - 5.48579909070e-4))
  
  
  expect_that(rdires$formula[rdorder] , equals(mtres$MF))
  
  
})

test_that("calcMF works with multiple mz inputs",{
  
  
  expect_that(calcMF(mz = 200.000659,
                     z = 1,
                     ppm = 5) , equals(calcMF(mz = c(200.000659,
                                                     200.000659),
                                              z = 1,
                                              ppm = 5)[[1]]))
  
  
})

test_that("calcMF works with parallel processing",{
  
  BiocParallel::register(BiocParallel::bpstart(
    if(Sys.info()['sysname'] == "Windows"){BiocParallel::SnowParam()
    }else{BiocParallel::MulticoreParam()}))
  
  
  expect_that(calcMF(mz = 200.000659,
                     z = 1,
                     ppm = 5) ,
              equals(calcMF(mz = c(200.000659,
                                   200.000659),
                            z = 1,
                            ppm = 5,
                            BPPARAM = BiocParallel::bpparam())[[1]]))
  
  
})


test_that("calcMF summarize works correctly",{
  
  
  expect_that(paste(calcMF(mz = 200.000659,
                           z = 1,
                           ppm = 5)$MF,"(+1)", sep = "", collapse = "|") , 
              equals(calcMF(mz = 200.000659,
                            z = 1,
                            ppm = 5,
                            summarize = TRUE)))
  
  
})


test_that("calcMF filters correctly",{
  
  
})