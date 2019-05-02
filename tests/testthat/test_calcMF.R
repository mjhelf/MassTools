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
  
 rdires <- Rdisop::decomposeMass(200.000659,
                       z = 1,
                       ppm = 5,
                       mzabs = 0,
                       elements = Rdisop::initializeCHNOPS())
  
 rdorder <- order(abs(rdires$exactmass - 200.000659 - 5.48579909070e-4))
 

 expect_that(rdires$formula[rdorder] , equals(mtres$MF))
 
  
})
