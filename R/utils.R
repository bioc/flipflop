## These functions are used in the main script to 
## - create the splicing graph 
## - calculate the cost of each node
## - basic check on the bins

# look if a given bin is striclty inside another given bin
#bin.inside <- function(inside.bin, outside.bin){
#   test <- FALSE
#   length.diff <- length(outside.bin) - length(inside.bin)
#   if(length.diff>=2){
#      test <- any( sapply(1:(length.diff-1), FUN=function(ii) identical( inside.bin , outside.bin[(ii+1):(ii+length(inside.bin))] ) ) )
#   }
#   return(test)
#}

# match a vector of indexes (bin) to a list of vectors of indexes (alllist)
# stop the loop at first match and return the position
quickmatch <- function(bin, alllist){
   test <- FALSE
   ii <- 0
   for(testbin in alllist){
      ii <- ii+1
      if(identical(testbin,bin)){
         test <- ii
         break
      }
   }
   return(test)
}

# given a particular suffix, find matching prefix among all given bins
# return the exon number that have to be added (the exon right after the prefix)
find.prefix <- function(suffix, bin, tophat.exons){
   test <- 0
   if( identical(bin[1:length(suffix)], suffix) & length(bin)>length(suffix) ){
      test <- bin[length(suffix)+1] 
   }
   return(as.integer(test))
}

# calculate the cost of each bin
# ie effective length
calculate.cost <- function(len.exons, bin, readlen){
   debut <- bin[1]
   fin <- bin[length(bin)]
   if(length(bin)>1){
      lint <- 0
      if(length(bin)>2){
         interieur <- bin[c(-1,-length(bin))]
         lint <- sum(len.exons[interieur])
      }
      cost <- max(0,min(len.exons[debut],readlen-lint-1) + min(len.exons[fin],readlen-lint-1) - readlen + lint + 1)
   }
   if(length(bin)==1){
      cost <- max(len.exons[debut] - readlen + 1,0)
   }
   return(cost)
}

# define if a bin can start or and a path in the graph
startstop <- function(len.exons, bin, readlen){
   logiquestart <- logiquestop <- FALSE
   if(length(bin)>1){
      teststart <- sum(len.exons[bin[-length(bin)]])
      teststop <- sum(len.exons[bin[-1]])
      if(teststart<readlen){ logiquestart <- TRUE }
      if(teststop<readlen){ logiquestop <- TRUE }
   }
   if(length(bin)==1){
      logiquestart <- TRUE
      logiquestop <- TRUE
   }
   return(c(logiquestart, logiquestop))
}

startstop_paired <- function(bin, list.byexons){
   logiquestart <- logiquestop <- FALSE
   if(length(bin)>1){
      teststart <- bin[1:(length(bin)-1)]
      teststop <- bin[2:length(bin)]
      testmatch.start <- quickmatch(teststart, list.byexons[[teststart[1]]])
      testmatch.stop <- quickmatch(teststop, list.byexons[[teststop[1]]])
      if(!testmatch.start){ logiquestart <- TRUE }
      if(!testmatch.stop){ logiquestop <- TRUE }
   }
   if(length(bin)==1){
      logiquestart <- TRUE
      logiquestop <- TRUE
   }
   return(c(logiquestart, logiquestop))
}

lastexon <- function(bin){
   return(bin[length(bin)])
}

firstexon <- function(bin){
   return(bin[1])
}

findtss <- function(possiblestart, exon.tss, tophat.exons, len.exons){
   test <- FALSE
   ind.possible <- firstexon(possiblestart)
   fin.tss <- tophat.exons[exon.tss,2]
   debut.possible <- tophat.exons[ind.possible,1]
   if(ind.possible %in% exon.tss) {
      test <- TRUE
   } else {
      tss.keep <- exon.tss[exon.tss<ind.possible]
      if(length(tss.keep)>0){
         tss.close <- tss.keep[which.min(ind.possible-tss.keep)]
         gap <- sum(tophat.exons[(tss.close+1):ind.possible,1] - tophat.exons[tss.close:(ind.possible-1),2])
         if(gap==0){
            test <- TRUE
         }
      }
   }
   return(test)
}

findpas <- function(possiblestop, exon.pas, tophat.exons, len.exons){
   test <- FALSE
   ind.possible <- lastexon(possiblestop)
   if(ind.possible %in% exon.pas){
      test <- TRUE
   }
   if(match(ind.possible, exon.pas, nomatch=0)==0){
      pas.keep <- exon.pas[exon.pas>ind.possible]
      if(length(pas.keep)>0){
         pas.close <- pas.keep[which.min(pas.keep-ind.possible)]
         gap <- sum(tophat.exons[(ind.possible+1):pas.close,1] - tophat.exons[ind.possible:(pas.close-1),2])
         if(gap==0){
            test <- TRUE
         }
      }
   }
   return(test)
}
