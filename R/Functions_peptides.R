#' @title PeptideMF
#'
#' convert peptide sequences into an MFobjects
#' 
#' @param seq one or more peptide sequences in a character vector 
#' @param AAs molecular formulas for each amino acid letter, defaults to proteinogenic amino acids
#' 
#' @return a list of MFobjects
#' 
#' @export
PeptideMF <- function(seq, AAs = c("A" = "C3H7NO2",  "R" = "C6H14N4O2",  "N" = "C4H8N2O3",
                                   "D" = "C4H7NO4",  "C" = "C3H7NO2S",  "Q" = "C5H10N2O3",
                                   "E" = "C5H9NO4",  "G" = "C2H5NO2",  "H" = "C6H9N3O2",
                                   "I" = "C6H13NO2",  "L" = "C6H13NO2",  "K" = "C6H14N2O2",
                                   "M" = "C5H11NO2S",  "F" = "C9H11NO2",  "P" = "C5H9NO2",
                                   "S" = "C3H7NO3",  "T" = "C4H9NO3",  "W" = "C11H12N2O2", 
                                   "Y" = "C9H11NO3",  "V" = "C5H11NO2")){
  seq <- gsub("[[:space:]]","",seq)
  
  pepAAs <- strsplit(seq, split = "")
  Ecounts <- lapply(pepAAs, function(p){makeMF(paste0(paste(AAs[p], collapse = ""),paste(rep("H-2O-1",length(p)-1), collapse = "") ))})
  
  return(Ecounts)
}

#' countAAs
#'
#' count amino acids in a sequence
#' 
#' @param seq one or more peptide sequences in a character vector 
#' 
#' @return a list of named numeric vectors
#' 
#' @export
countAAs <- function(seq){
  
  seq <- gsub("[[:space:]]","",seq)
  
  pepAAs <- strsplit(seq, split = "")
  
  uniques <- lapply(pepAAs, unique)
  
  res <- mapply(function(u,p){
    counts <- sapply(u,function(aa){sum(p==aa)});
    names(counts) <- u;
    return(counts)},
    uniques,pepAAs, SIMPLIFY = F)
  
  return(res)
}


#' permutatePeptideMass
#'
#' Apply a set of mass modifications to a peptide with a known mass
#'
#' @param df data.frame with at least two columns, one with a non-empty 
#' sequence string and one with a mass for that sequence
#' @param modifications a data.frame with colums \code{AAs} 
#' @param sequenceCol name of the sequence column in df
#' @param massCol name of the mass column in df
#' @param chargeCol name of column that specifies charge state (used to adjust 
#' modification masses)
#' 
#' @details 
#' \subsection{modifications}{
#' Columns that need to be in the \code{modifications data.frame}:
#' \describe{
#' \item{AAs}{containing letters for the amino acids (e.g. "TS")
#' affected by this modification}
#' \item{min}{minimum number of times this modification can be applied}
#' \item{max}{maximum number of times this modification can be applied}
#' \item{mass}{mass of this modification}
#' \item{tag}{ label for this modification (will automatically be preceded by 
#' number of these modifications applied)}
#' }
#' min and max can be positive, negative or 0. The larger *absolute* value has 
#' to be in the \code{max} column
#' }
#' 
#' @examples
#' peptides = data.frame(seq = "ASDFGHKLIPYTREWQCNM",mz = 2295.056)
#' mods = data.frame(AAs = c("TS", ""),
#'  min = c(0,1), max = c(-10,2),
#'   mass = c(18.01,12.0), tag = c("DH","C") )
#' 
#' permutatePeptideMass(peptides, mods)
#'
#' 
#' @seealso \code{\link{permutateMass}}
#' 
#' @export
permutatePeptideMass <- function(df, modifications,
                                 sequenceCol = "seq", 
                                 massCol = "mz",
                                 chargeCol = NULL){
  if((missing(modifications) 
      || is.null(modifications) 
      ||!nrow(modifications))
     ){
    
    if(!"modifications" %in% colnames(df)){
   df$modifications <- ""
    }
   
   return(df)
   
  }
  
  res <-  lapply(seq(nrow(df)),function(i, modifications){
    
    AAcounts <- countAAs(df[[sequenceCol]][i])[[1]]
    
    selective <- which(modifications$AAs != "" & !is.na(modifications$AAs))
    
    if(length(selective) > 0){
      
      selectivecounts <- sapply(selective,function(sel){
        
        targets <- countAAs(modifications$AAs[sel])[[1]]
        
        targetAAs <- match(names(targets), names(AAcounts))
        
        if(!any(!is.na(targetAAs))){
          return(0)
        }
        return(sum(AAcounts[na.omit(targetAAs)]))
        
      })
      
      #make sure the minimum and maximum values are in line with the number of required amino acids, and allow the numbers to be negative,as well
      modifications$min[selective] <- sapply(seq(length(selective)),
                                             function(i){ sign(modifications$min[selective][i]) * min(c(abs(modifications$min[selective][i]),
                                                                                                        selectivecounts[i])) })
      
      modifications$max[selective] <- sapply(seq(length(selective)),
                                             function(i){ sign(modifications$max[selective][i]) * min(c(abs(modifications$max[selective][i]),
                                                                                                        selectivecounts[i])) })
      
    }
    
    #this means that users HAVE to put the larger *absolute* value to into the max column
    #would advise against mixing negative and positive values in min and max columns!
    modifications <- modifications[abs(modifications$max) >= abs(modifications$min),]
    
    #use charge state
    if(!is.null(chargeCol)){
      modifications$mass <- modifications$mass/df[[chargeCol]][i]
    }
    permutations <- permutateMass(df[[massCol]][i],modifications)
    
    colnames(permutations)[colnames(permutations) == "mass"] <- massCol
    
    suppressWarnings({
      merged <- cbind(df[i, colnames(df) != massCol, drop = F], permutations)
    })
    
    return(merged)
    
  }, modifications = modifications)
  
  return(do.call(rbind,res))
}


#' plotAnnotatedPeptide
#'
#' plot an annotated peptide sequence
#'
#' @param peaklabels data.frame with labeling informationi,
#'  as returned by \code{\link{annotateSpectrum}}
#' @param sequence peptide sequence (character)
#' @param parseLabels optionally apply parse() to  
#' peaklabels$label and add as label (NYI)
#' @param sortby sort (decreasing) by this column. Only the top 2 ions 
#' for each position (for both a,b,c ions AND x,y,z ions) will be shown
#' @param cx font size factor for the sequence characters,
#'  other font sizes will be scaled accordingly
#' @param yoffset y axis offset
#'
#' @export
plotAnnotatedPeptide <- function(peaklabels,
                                 sequence,
                                 parseLabels = F,
                                 sortby = "intensitySpec",
                                 cx = 2,
                                 yoffset = 0){
  
  
  
  par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i", xpd = NA)
  
  chars <- strsplit(sequence, split = "")[[1]]
  
  plot(numeric(),
       numeric(),
       ylim = c(-1,1),
       xlim = c(0,1),
       xaxs = "i", yaxs = "i",
       type = "n", ann = FALSE, bty = "n",
       axes = F#, asp =0.5
  )
  
  charheight <- strheight("M", cex = cx)
  charwidth <- strwidth("M", cex = cx)
  
  ionheight <- strheight("M", cex = 0.5*cx)
  ionwidth <- strwidth("M", cex = 0.5*cx)
  
  #this is now used as default, e.g. by strwidth() calls
  par(cex = 0.5*cx)
  
  spacer <- 0.45*ionheight
  bigspacer <- 2*spacer
  
  #conversionfactor between strheight and strwidth, assuming M is a square 
  #(tested M, and it is a square in terms of plotting dimensions)
  HtoWratio <- charheight/charwidth
  
  fieldwidth <- charwidth*1.02
  centers <- (seq(length(chars))*fieldwidth-0.5*fieldwidth  
              + (0.5-0.5*fieldwidth*length(chars)))
  #last part centers sequence in plot area
  
  bottoms <- -0.55*charheight + yoffset
  tops <- 0.5*charheight + yoffset
  bionpos <-  tops + 0.5*charheight
  yionpos <- bottoms - 0.5*charheight
  
  
  text(centers,yoffset,labels=chars, cex=cx, adj=c(0.5,0.5))
  
  if(nrow(peaklabels)>0){
    
    if(sortby %in% colnames(peaklabels)){
      
      peaklabels <- peaklabels[order(peaklabels[[sortby]],
                                     decreasing = T),]
      
    }
    
    peaklabels$ppm <- ((peaklabels$mz-peaklabels$mzSpec)/peaklabels$mzSpec)*1e6
    
    massionlabels <-  format(round(peaklabels$mzSpec,5),
                             nsmall = 5, scientific = F)
    
    massionlabels[peaklabels$z == 1] <- paste0(massionlabels[peaklabels$z == 1],
                                               "^{+{}}")
    
    massionlabels[peaklabels$z > 1] <- paste0(massionlabels[peaklabels$z > 1],
                                              "^{",
                                              peaklabels$z[peaklabels$z > 1],
                                              "+{}}")
    
    
    
    #adjust the y ion positions to same scale as b ions
    peaklabels$pos[peaklabels$type %in% c("x","y","z")] <- length(chars) - peaklabels$pos[peaklabels$type %in% c("x","y","z")]
    
    #for drawing the fragmentation lines
    presentB <- unique(peaklabels$pos[peaklabels$type %in% c("a","b","c")])
    presentY <- unique(peaklabels$pos[peaklabels$type %in% c("x","y","z")])
    alllines <- unique(c(presentB,presentY))    
    
    
    
    
    #vertical lines
    segments(centers[alllines]+0.5*fieldwidth,rep(bottoms,
                                                  length(alllines)),
             centers[alllines]+0.5*fieldwidth,rep(tops, 
                                                  length(alllines)),
             col = "black", lty = "solid", lwd = 1)
    
    segments(centers[presentB]+0.5*fieldwidth,rep(tops, 
                                                  length(presentB)),
             centers[presentB]+0.1*fieldwidth,rep(bionpos - 0.05*charheight,
                                                  length(presentB)),
             col = "black", lty = "solid", lwd = 1)
    
    segments(centers[presentY]+0.5*fieldwidth,rep(bottoms,
                                                  length(presentY)),
             centers[presentY]+0.9*fieldwidth,rep(yionpos + 0.05*charheight,
                                                  length(presentY)),
             col = "black", lty = "solid", lwd = 1)
    
    #this way of selection treats only picks the first one AND treats 
    #a b c ions the same (as opposed to just !duplicated(peaklabels$ion) 
    bsel <- which(peaklabels$type %in% c("a","b","c"))
    bsel <- bsel[!duplicated(peaklabels$pos[bsel])]
    
    ysel <- which(peaklabels$type %in% c("x","y","z"))
    ysel <- ysel[!duplicated(peaklabels$pos[ysel])]
    
    if(parseLabels){
      labelsfi <- parse(text = peaklabels$label)
    }else{
      labelsfi <- peaklabels$mz
    }
    
    parsedlabs <- parse(text = gsub("\\[\\]","",paste0('"[',
                                                       peaklabels$modifications,
                                                       ']"')))
    
    #now find the second ranking ions
    if(length(bsel)>0){
      bsel2 <- which(peaklabels$type %in% c("a","b","c"))
      bsel2 <- bsel2[!bsel2 %in% bsel]
      bsel2 <- bsel2[!duplicated(peaklabels$pos[bsel2])]
      
      
      #label the b ions
      text(centers[presentB],
           rep(bionpos, length(bsel)),
           labels = parse(text = paste0("b[",presentB,"]^{}")),
           adj=c(0.5,0), srt=0, cex = 0.5*cx)
      
      
      
      ##First rank labels
      
      ### b ion modifications
      
      
      
      bmodlabels <- parsedlabs[bsel]
      bion1modpos <- bionpos + 1.55*ionheight
      
      text(centers[peaklabels$pos[bsel]],
           rep(bion1modpos, length(bsel)),
           labels = bmodlabels,
           adj=c(1,0.5), srt=270, cex = 0.5*cx)
      
      #ppm labels
      bppmpos <-  (bion1modpos
                   + bigspacer 
                   + max(strwidth(bmodlabels))*HtoWratio)
      bppmlabels <- paste0("(",format(round(peaklabels$ppm,1)[bsel],
                                      nsmall = 1, scientific = F),")")
      
      text(centers[peaklabels$pos[bsel]],
           rep(bppmpos, length(bsel)),
           labels = bppmlabels,
           adj=c(1,0.5), srt=270, cex = 0.5*cx)
      
      
      
      bmasspos <-  (bppmpos 
                    + spacer 
                    + max(strwidth(bppmlabels))*HtoWratio)
      bmasslabels <- parse(text = massionlabels[bsel])
      
      text(centers[peaklabels$pos[bsel]],
           rep(bmasspos, length(bsel)),
           labels = bmasslabels,
           adj=c(1,0.5), srt=270, cex = 0.5*cx)
      
      
      
      ##Second rank labels
      
      ### b ion modifications
      
      if(length(bsel2)){
        
        bmodlabels <- parsedlabs[bsel2]
        bion1modpos <- (bmasspos 
                        + bigspacer 
                        + max(strwidth(bmasslabels))*HtoWratio)
        
        text(centers[peaklabels$pos[bsel2]],
             rep(bion1modpos, length(bsel2)),
             labels = bmodlabels,
             adj=c(1,0.5), srt=270, cex = 0.5*cx)
        
        bppmpos <-  (bion1modpos
                     + bigspacer
                     + max(strwidth(bmodlabels))*HtoWratio) 
        bppmlabels <- paste0("(",format(round(peaklabels$ppm,1)[bsel2],
                                        nsmall = 1, scientific = F),")")
        
        text(centers[peaklabels$pos[bsel2]],
             rep(bppmpos, length(bsel2)),
             labels = bppmlabels,
             adj=c(1,0.5), srt=270, cex = 0.5*cx)
        
        bmasspos <-  (bppmpos 
                      + spacer 
                      + max(strwidth(bppmlabels))*HtoWratio) 
        bmasslabels <- parse(text = massionlabels[bsel2])
        
        text(centers[peaklabels$pos[bsel2]],
             rep(bmasspos, length(bsel2)),
             labels = bmasslabels,
             adj=c(1,0.5), srt=270, cex = 0.5*cx)
        
        
        
        
        
      }
      
    }
    
    if(length(ysel)){
      ysel2 <- which(peaklabels$type %in% c("x","y","z"))
      ysel2 <- ysel2[!ysel2 %in% ysel]
      ysel2 <- ysel2[!duplicated(peaklabels$pos[ysel2])]
      
      #label the y ions
      text(centers[presentY]+1*fieldwidth,
           rep(yionpos, length(ysel)),
           labels = parse(text = paste0("y[",
                                        length(chars)-presentY,
                                        "]^{}")),
           adj=c(0.5,1), srt=0, cex = 0.5*cx)
      
      ### y ion modifications
      ymodlabels <- parsedlabs[ysel]
      yion1modpos <- yionpos - 1.5*ionheight
      
      text(centers[peaklabels$pos[ysel]]+fieldwidth,
           rep(yion1modpos, length(ysel)),
           labels = ymodlabels,
           adj=c(1,0.5), srt=90, cex = 0.5*cx)
      
      yppmpos <-  (yion1modpos 
                   - bigspacer 
                   - max(strwidth(ymodlabels, cex = 0.5*cx))*HtoWratio)
      yppmlabels <-paste0("(",format(round(peaklabels$ppm,1)[ysel],
                                     nsmall = 1, scientific = F),")")
      
      text(centers[peaklabels$pos[ysel]]+fieldwidth,
           rep(yppmpos, length(ysel)),
           labels = yppmlabels,
           adj=c(1,0.5), srt=90, cex = 0.5*cx)
      
      ymasspos <-  (yppmpos 
                    - spacer 
                    - max(strwidth(yppmlabels))*HtoWratio)
      ymasslabels <- parse(text = massionlabels[ysel])
      
      text(centers[peaklabels$pos[ysel]]+fieldwidth,
           rep(ymasspos, length(ysel)),
           labels = ymasslabels,
           adj=c(1,0.5), srt=90, cex = 0.5*cx)
      
      if(length(ysel2)){
        
        ymodlabels <-parsedlabs[ysel2]
        yion1modpos <- (ymasspos 
                        - bigspacer 
                        - max(strwidth(ymasslabels))*HtoWratio)
        
        text(centers[peaklabels$pos[ysel2]]+fieldwidth,
             rep(yion1modpos, length(ysel2)),
             labels = ymodlabels,
             adj=c(1,0.5), srt=90, cex = 0.5*cx)
        
        yppmpos <-  (yion1modpos
                     - bigspacer 
                     - max(strwidth(ymodlabels))*HtoWratio)
        yppmlabels <-paste0("(",format(round(peaklabels$ppm,1)[ysel2],
                                       nsmall = 1, scientific = F),")")
        
        text(centers[peaklabels$pos[ysel2]]+fieldwidth,
             rep(yppmpos, length(ysel2)),
             labels = yppmlabels,
             adj=c(1,0.5), srt=90, cex = 0.5*cx)
        
        ymasspos <-  (yppmpos
                      - spacer
                      - max(strwidth(yppmlabels))*HtoWratio)
        ymasslabels <- parse(text = massionlabels[ysel2])
        
        text(centers[peaklabels$pos[ysel2]]+fieldwidth,
             rep(ymasspos, length(ysel2)),
             labels = ymasslabels,
             adj=c(1,0.5), srt=90, cex = 0.5*cx)
        
        
        
        
        
      }
      
    }
  }
  
  #make sure the custom cex in function doesn't interfere
  par(cex = 1)
  
}