## These functions are used in the main script to 
## - create the splicing graph 
## - calculate the cost of each node
lastexon <- function(namebin){
   attrape <- as.numeric(strsplit(namebin, split='-')[[1]])
   ind.dernier <- length(attrape)
   return(attrape[ind.dernier])
}
firstexon <- function(namebin){
   attrape <- as.numeric(strsplit(namebin, split='-')[[1]])
   return(attrape[1])
}
find.prefix <- function(suffix, datatype, tophat.exons){
   tocompare <- rep(0,suffix[length(suffix)])
   tocompare[suffix] <- 1
   test <- (identical(datatype[1:suffix[length(suffix)]], tocompare) &  sum(datatype[(suffix[length(suffix)]+1):length(datatype)])>=1)
   dernier <- suffix[length(suffix)]
   premier <- min(which(datatype==1))
   if((dernier+1)==premier & tophat.exons[dernier,2]==tophat.exons[premier,1]){ # test is true for neighbours exons
      test <- TRUE
   }
   return(test)
}
calculate.cost <- function(len.exons, namebin, readlen){
   attrape <- as.numeric(strsplit(namebin, split='-')[[1]])
   debut <- attrape[1]
   fin <- attrape[length(attrape)]
   if(length(attrape)>1){
      lint <- 0
      if(length(attrape)>2){
         interieur <- attrape[c(-1,-length(attrape))]
         lint <- sum(len.exons[interieur])
      }
      cost <- max(0,min(len.exons[debut],readlen-lint-1) + min(len.exons[fin],readlen-lint-1) - readlen + lint + 1)
   }
   if(length(attrape)==1){
      cost <- max(len.exons[debut] - readlen + 1,0)
   }
   return(cost)
}
startstop <- function(len.exons, namebin, readlen, count){
   attrape <- as.numeric(strsplit(namebin, split='-')[[1]])
   logiquestart <- logiquestop <- FALSE
#   if(length(attrape)>1 & count[namebin,]!=0){
   if(length(attrape)>1){
      teststart <- sum(len.exons[attrape[-length(attrape)]])
      teststop <- sum(len.exons[attrape[-1]])
      if(teststart<readlen){ logiquestart <- TRUE }
      if(teststop<readlen){ logiquestop <- TRUE }
   }
#   if(length(attrape)==1 & count[namebin,]!=0){
   if(length(attrape)==1){
      logiquestart <- TRUE
      logiquestop <- TRUE
   }
   return(c(logiquestart, logiquestop))
}
startstop_paired <- function(namebin, allbins){
   attrape <- as.numeric(strsplit(namebin, split='-')[[1]])
   logiquestart <- logiquestop <- FALSE
   if(length(attrape)>1){
      teststart <- paste(attrape[1:(length(attrape)-1)], collapse='-')
      teststop <- paste(attrape[2:length(attrape)], collapse='-')
      if(match(teststart, allbins, nomatch=0)==0){ logiquestart <- TRUE }
      if(match(teststop, allbins, nomatch=0)==0){ logiquestop <- TRUE }
   }
   if(length(attrape)==1){
      logiquestart <- TRUE
      logiquestop <- TRUE
   }
   return(c(logiquestart, logiquestop))
}
findtss <- function(possiblestart, tss.find, tophat.exons, len.exons){
   test <- FALSE
   ind.possible <- firstexon(possiblestart)
   fin.tss <- tophat.exons[tss.find,2]
   debut.possible <- tophat.exons[ind.possible,1]
   if(ind.possible %in% tss.find){
      test <- TRUE
   }
   if(match(ind.possible, tss.find, nomatch=0)==0){
      tss.keep <- tss.find[tss.find<ind.possible]
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
findpas <- function(possiblestop, pas.find, tophat.exons, len.exons){
   test <- FALSE
   ind.possible <- lastexon(possiblestop)
   if(ind.possible %in% pas.find){
      test <- TRUE
   }
   if(match(ind.possible, pas.find, nomatch=0)==0){
      pas.keep <- pas.find[pas.find>ind.possible]
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

