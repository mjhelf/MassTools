context("Plotting functions")


test_that("peptide annotation and plotting works",{

  bm <- 2295.056
  peptides = data.frame(seq = "ASDFGHKLIPYTREWQCNM",mz = bm, stringsAsFactors = F)
  mods = data.frame(AAs = c("TS", ""), min = c(0,1), max = c(-10,2), mass = c(18.01,12.0), tag = c("DH","C") )
 
  # Use these lines to remake the base fragments. using .rds file to avoid dependency for MSnbase (complicates Travis build)
  # unmodfrags <- MSnbase::calculateFragments(sequence = peptides$seq[1],
  #                                           neutralLoss = list(water = character(0), ammonia = character(0)),
  #                                           z = 1)
  #   saveRDS(unmodfrags, "testFragments.rds")

    unmodfrags <- readRDS("testFragments.rds")
  
  fragments <- permutatePeptideMass(unmodfrags, mods )
  
  annotation <- annotateSpectrum(fragments, data.frame(mz = fragments$mz,
                                                       intensity = 1000+  seq_along(fragments$mz)),
                                             mzlabel = F,
                                             unmodifiedTag = "",
                                             ppm = 5,
                                             abs = 0)
  
  expect_equivalent(annotation[1:3,],
                    data.frame(mzSpec = c(84.044386, 96.044386, 171.076416),
                               intensitySpec= c(1001, 1002, 1003),
                               ion = c("b1","b1","b2"),
                               type = c("b","b","b"),
                               pos = c(1,1,2),
                               z = 1,
                               seq = c("A","A","AS"),
                               mz = c(84.044386, 96.044386, 171.076416),
                               modifications = c('1C','2C','1C'),
                               color = "blue3",
                               label = c('b[1]^{+{}}* "[1C]"','b[1]^{+{}}* "[2C]"', 'b[2]^{+{}}* "[1C]"'),
                               stringsAsFactors = FALSE
                               ))
  
  
  #testing different ways of having no modifications:

  expect_equal(permutatePeptideMass(unmodfrags, modifications = mods[logical(0),] ),
               permutatePeptideMass(unmodfrags, modifications = NULL))
  
  expect_equal(permutatePeptideMass(unmodfrags, modifications = mods[logical(0),] ),
               permutatePeptideMass(unmodfrags))
  
  unmodfrags2 <- unmodfrags
  unmodfrags2$modifications <- ""
  
  expect_equal(permutatePeptideMass(unmodfrags, modifications = mods[logical(0),] ),
               unmodfrags2)
  

  annotation2 <- annotateSpectrum(unmodfrags2, data.frame(mz = unmodfrags2$mz,
                                                       intensity = 1000+  seq_along(unmodfrags2$mz)),
                                 mzlabel = F,
                                 unmodifiedTag = "",
                                 ppm = 5,
                                 abs = 0)
  
  
  expect_doppelganger("no_modifications",
                      plotAnnotatedPeptide(peaklabels = annotation2,
                                           sequence = "ASDFGHKLIPYTREWQCNM",
                                           parseLabels = T,
                                           sortby = "intensitySpec",
                                           cx = 2,
                                           yoffset = 0))
  
  expect_doppelganger("with_modifications",
                      plotAnnotatedPeptide(peaklabels = annotation,
                                           sequence = "ASDFGHKLIPYTREWQCNM",
                                           parseLabels = T,
                                           sortby = "intensitySpec",
                                           cx = 2,
                                           yoffset = 0))

})
