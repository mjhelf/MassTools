context("Spectra comparison functions")


test_that("pairCompare works",{
  
  expect_equal(MassTools:::pairCompare(c(1,0,1),c(1,0,1),
                          NAasZero = T,
                          method = c("cosine")),1)
  
  expect_equal(MassTools:::pairCompare(c(1,0,1),c(0,1,0),
                          NAasZero = T,
                          method = c("cosine")),0)
  
  expect_equal(MassTools:::pairCompare(c(1,1),c(1,1),
                                       NAasZero = T,
                                       method = c("cosine")),1)
  
  expect_equal(MassTools:::pairCompare(NA,NA,
                                       NAasZero = T,
                                       method = c("cosine")),NaN)
  
  expect_equal(MassTools:::pairCompare(1,1,
                                       NAasZero = T,
                                       method = c("cosine")),1)
  
  expect_equal(MassTools:::pairCompare(c(1,NA,1),c(1,NA,1),
                                       NAasZero = T,
                                       method = c("cosine")),1)
  
  expect_equal(MassTools:::pairCompare(c(1,0,1),c(1,NA,1),
                                       NAasZero = T,
                                       method = c("cosine")),1)
  
  
  
  #pearson method is handling some edge cases in suboptimal ways 
  #for large scale usage
  expect_equal(MassTools:::pairCompare(c(1,0,1),c(1,0,1),
                          NAasZero = T,
                          method = c("pearson")),1)
  
  
  expect_equal(MassTools:::pairCompare(c(1,0,1),c(0,1,0),
                          NAasZero = T,
                          method = c("pearson")),-1)
  
  expect_warning(expect_equal(MassTools:::pairCompare(c(1,1),c(1,1),
                                       NAasZero = T,
                                       method = c("pearson")),NA_real_))
  
 expect_equal(MassTools:::pairCompare(1,1,
                                       NAasZero = T,
                                       method = c("pearson")),NA_real_)
 
  expect_warning(expect_equal(MassTools:::pairCompare(c(1,NA,1),c(1,NA,1),
                          NAasZero = T,
                          method = c("pearson")), NA_real_))
  
  
})

test_that("network1 and makeEdges work",{

    spec1 <- matrix(c(100, 115, 150,
                    1000, 2000, 1000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
    
    spec2 <- matrix(c(100, 115, 150,
                    1000, 500, 1000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
    
    spec3 <- matrix(c(115, 1501,
                      500, 1000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
    
    spec4 <- matrix(c(115, 150,
                      500, 1000), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
    
    multispec <- lapply(1:10,function(x){
        matrix(c(100,  150,
                 1000, 1000*x), ncol = 2, dimnames = list(NULL, c("mz","intensity")))
        
        })
    
    expect_equal(network1(spec1, spec1, mztol = 0.005,
                         parentshift = 0,
                         method = "cosine",
                         minpeaks = 3,
                         nonmatched = T),1)
    
    expect_equal(network1(spec1, spec1[FALSE,], mztol = 0.005,
                          parentshift = 0,
                          method = "cosine",
                          minpeaks = 3,
                          nonmatched = T),NA_real_)
    
    expect_equal(network1(spec1, spec1+0.005, mztol = 0.005,
                          parentshift = 0,
                          method = "cosine",
                          minpeaks = 3,
                          nonmatched = T),1)
    
    expect_equal(network1(spec1, spec1+0.006, mztol = 0.005,
                          parentshift = 0,
                          method = "cosine",
                          minpeaks = 3,
                          nonmatched = T),0)
    
    expect_equal(network1(spec1, spec1+0.006, mztol = 0.005,
                          parentshift = 0.006,
                          method = "cosine",
                          minpeaks = 3,
                          nonmatched = T),1)
    
    #testing behavior if one or more peaks from one spectrum match one or more from the other...
    expect_equal(network1(spec1, rbind(spec1,spec1), mztol = 0.05,
                                     parentshift = 0,
                                     method = "cosine",
                                     minpeaks = 4,
                                     nonmatched = T),0)
    
    expect_equal(network1(rbind(spec1,spec1), rbind(spec1,spec1), mztol = 0.05,
                          parentshift = 0,
                          method = "cosine",
                          minpeaks = 4,
                          nonmatched = T),0)
    
    expect_equal(network1(rbind(spec1,spec1), rbind(spec1,spec1), mztol = 0.05,
                          parentshift = 0,
                          method = "cosine",
                          minpeaks = 3,
                          nonmatched = T),1)
    
    expect_equal(network1(rbind(spec1,spec1), spec1, mztol = 0.05,
                          parentshift = 0,
                          method = "cosine",
                          minpeaks = 3,
                          nonmatched = T),1)
    expect_equal(network1(rbind(spec1,spec1), spec1, mztol = 0.05,
                          parentshift = 0,
                          method = "cosine",
                          minpeaks = 4,
                          nonmatched = T),0)
    
    expect_equal(network1(spec1, rbind(spec1+c(1:3,0,0,0),spec1), mztol = 2,
                                     parentshift = 0,
                                     method = "cosine",
                                     minpeaks = 3,
                                     nonmatched = F),1)
    
    
    
    expect_equal(round(network1(spec1, spec2, mztol = 0.005,
             parentshift = 0,
             method = "cosine",
             minpeaks = 3,
             nonmatched = T),3),0.816)
    
    expect_equal(round(network1(spec1, spec4, mztol = 0.005,
                                parentshift = 0,
                                method = "cosine",
                                minpeaks = 3,
                                nonmatched = T),3),0)
    
    
    #effect of nonmatched argument:
    expect_equal(round(network1(spec1, spec4, mztol = 0.005,
                                parentshift = 0,
                                method = "cosine",
                                minpeaks = 2,
                                nonmatched = T),3),0.730)
    
    expect_equal(round(network1(spec2, spec4, mztol = 0.005,
                                parentshift = 0,
                                method = "cosine",
                                minpeaks = 2,
                                nonmatched = F),3),1)
    
    
    expect_equal(makeEdges(list(spec1, spec1, spec2),
              parentmasses = NULL,
              mztol = 0.005, 
              method = "cosine",
              minpeaks = 3, 
              nonmatched = T),data.frame(from = c(1,1,2),
                                         to = c(2,3,3),
                                         cosine = c(1,0.816496576641330,0.816496576641330)))
    
    expect_equal(makeEdges(list(spec1, spec1, spec2),
                           parentmasses = c(1,1,2),
                           mztol = 0.005, 
                           method = "cosine",
                           minpeaks = 3, 
                           nonmatched = T),data.frame(from = c(1,1,2),
                                                      to = c(2,3,3),
                                                      cosine = c(1,0.816496576641330,0.816496576641330),
                                                      deltamz = c(0,-1,-1)))
        
    })
    