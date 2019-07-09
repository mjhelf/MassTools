context("Calculating mz values")

test_that("calcIons correctly calculates ion m/z values",{
  
  ions <- calcIons(c("CH4", "H2O"),charges = c(-1, 1, 0, -2, 2))
  
  chemcalc_results <- c("CH4-1" = 15.02402,
                        "H2O-1" = 17.00329,
                        "CH4+1" = 17.03858,
                        "H2O+1" = 19.01784,
                        "CH4" = 16.03130, 
                        "H2O" = 18.01056,
                        "CH4-2" = 7.00837,
                        "H2O-2" = 7.99801,
                        "CH4+2" = 9.02293,
                        "H2O+2" = 10.01256)
  
  ref_df <- data.frame(mz = chemcalc_results,
                       formula = rep(c("CH4","H2O"),5),
                       charge = c(-1,-1,1,1,0,0,-2,-2,2,2),
                       ion = c("[M-H]1-", "[M-H]1-",
                               "[M+H]1+", "[M+H]1+",
                               "[M]", "[M]",
                               "[M-2H]2-", "[M-2H]2-",
                               "[M+2H]2+", "[M+2H]2+"),
                       stringsAsFactors = F)
  
  ions$mz <- round(ions$mz,5)
  
  expect_equivalent(ions, ref_df)
  
  ions2 <- calcIons(c("CH4NOPSNa"),charges = c(-1, 0, 2))
  ions2$mz <- round(ions2$mz,5)
  
  chemcalc_results2 <- c("CH4NOPSNa-1" = 130.95761,
                         "CH4NOPSNa" = 131.96489,
                         "CH4NOPSNa+2" = 66.98972
                         )
  
  ref_df2 <- data.frame(mz = chemcalc_results2,
                       formula = rep(c("CH4NOPSNa"),3),
                       charge = c(-1,0,2),
                       ion = c("[M-H]1-",
                               "[M]",
                               "[M+2H]2+"),
                       stringsAsFactors = F)
  
  expect_equivalent(ions2, ref_df2)
  
  
  
  })

