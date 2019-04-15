context("Building MFobjects")

test_that("makeMF works correctly",{
  
  water <- getOption("MassTools.elements")
  
  water["H"] <- 2
  water["O"] <- 1
  
 expect_that(unclass(makeMF("H2O")), equals(water))
 expect_that(unclass(makeMF("")), equals(getOption("MassTools.elements")))
 expect_that(unclass(makeMF(paste0(paste(names(water), collapse = "")))), equals(getOption("MassTools.elements")+1))
 expect_that(unclass(makeMF(paste0(paste(names(water), collapse = "2"),"2"))), equals(getOption("MassTools.elements")+2))
 expect_that(unclass(makeMF(paste0(paste(names(water),names(water), collapse = "", sep = "")))), equals(getOption("MassTools.elements")+2))
 
 #MF string with spaces:
 expect_that(unclass(makeMF(paste0(paste(names(water),names(water), collapse = "", sep = " ")))), equals(getOption("MassTools.elements")+2))
 
  
})