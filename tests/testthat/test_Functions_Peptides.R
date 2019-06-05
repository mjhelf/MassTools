context("Peptide functions")

test_that("Molecular formulas of peptides are calculated correctly",{
  
 expect_match(remakeMF(PeptideMF("ASDFGHKLIPYTREWQCNM")[[1]]),
              "C102H150N28O29S2")
  
  #test length >1
  expect_true(compare(sapply(PeptideMF(c("ASDFGHKLIPYTREWQCNM",
                                         "ASDFGHKLIPYTREWQCNM")), remakeMF),
                 c("C102H150N28O29S2","C102H150N28O29S2"))$equal)
  
  #test whitespace removal
  expect_match(remakeMF(PeptideMF("ASDFGHK
                                  LIPY  TREWQC  NM")[[1]]),
               "C102H150N28O29S2")
  
})





test_that("counting AAs in peptide sequence works",{
  
  
  expect_equal(countAAs("AGAGA")[[1]],
               c(A=3,G=2))
  
  #test length >1
  expect_equal(countAAs(c("AGAGA",
                          "SDFGSDFI")),
               list(c(A=3,G=2),
                    c(S=2,D=2,F=2,G=1,I=1)))
  
  #test whitespace removal
  expect_equal(countAAs("AG
                        A   G  A")[[1]],
               c(A=3,G=2))
  
  })

test_that("mergeMS counts combined peaks",{

  bm <- 2295.056
  peptides = data.frame(seq = "ASDFGHKLIPYTREWQCNM",mz = bm)
  mods = data.frame(AAs = c("TS", ""), min = c(0,1), max = c(-10,2), mass = c(18.01,12.0), tag = c("DH","C") )
 
  expect_equal(permutatePeptideMass(peptides, mods, 
                                    massCol = "mz",
                                    sequenceCol = "seq"
                                    )$mz,
               c(bm+12, bm+24,bm+12-18.01, bm+24-18.01,bm+12-36.02, bm+24-36.02)
  )
  
  
})
