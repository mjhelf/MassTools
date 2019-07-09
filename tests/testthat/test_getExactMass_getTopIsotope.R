context("Calculating masses with getExactMass and getTopIsotope")

test_that("getExactMass correctly calculates exact masses",{
  
  
  expect_equivalent(getExactMass(c("CH4", "C100H200", "C10000H20000")),
                    c(16.0313001283, 1401.5650064142, 140156.50064142))
  
  expect_equivalent(getExactMass(c("C-1H-2", "PO4NaH-2")),
               c(-14.0156500641, 117.9431893912 - 2.0156500641))
  
})

test_that("getTopIsotope correctly provides mass of most abundant isotope",{
  
  #note: first two reference isotopes calculated with ecipex package
  expect_equivalent(getTopIsotope(c("CH4", "C100H200", "C1000H2000")),
                    c(16.0313001283, 1402.5683612518, 14026.68782535))
  
  #NOTE: For "C1000H2000", both enviPat and ecipex calculate 14025.68361 as most abundant,
  # as opposed to getTopIsotope's 14026.68782535
  
  expect_equivalent(getTopIsotope(c("C-1H-2", "PO4NaH-2")),
                    c(-14.0156500641, 117.9431893912 - 2.0156500641))
  
})
