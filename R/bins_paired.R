#------ This function elaborate the bins when reads are paired-end -----# 
# Bins are the node of the splicing graph. 


bins_paired <- function(inpf, npetype, tophat.exons, count.exons, len.exons, n.exons, readtype, gap, std, n.samples.max, samples.name, n.samples) {
   
   longtype <- vector(mode='list') # record each created type
   numbertype <- vector(mode='list') # record the number of each type
   compteur <- 0

   numbertype.samples <- vector(mode='list') # record the number of each type for each sample
   count.first.samples <- NULL

   # I) first pass - create types based on left/right mapping
   while(compteur<npetype){
      compteur <- compteur+1
      plouf <- scan(inpf, what=character(0), nlines=2, quiet=TRUE)
      sgone <- as.numeric(plouf[1]) # left type index
      sgtwo <- as.numeric(plouf[2]) # right type index
      nimp <- as.numeric(plouf[3]) # number of such types
     
      oneleft <- which(readtype[sgone,1:n.exons]==1) # type of left read
      oneright <- which(readtype[sgtwo,1:n.exons]==1) # type of right read
      lastleft <- oneleft[length(oneleft)]
      firstright <- oneright[1]
      
      implicated <- as.integer(unlist(strsplit(plouf[4:(4+as.numeric(plouf[3])-1)],split=':')))
      npair <- implicated[seq(2,2*nimp,2)] # number of such paired-end reads
      alldist <- implicated[seq(1,2*nimp,2)] # all insert distances
     
      # read the lines specific to each sample
      if(n.samples.max > 1) {
         plouf.samples <- scan(inpf, what=character(0), nlines=n.samples.max, quiet=TRUE)
      }
      if(!is.null(samples.name)){
         info.samples <- lapply(samples.name, FUN=function(sn){
                                aa <- which(plouf.samples==sn)
                                nimp <- as.integer(plouf.samples[aa+1])
                                ns <- 0
                                ds <- NULL
                                if(nimp>0){
                                   imp <- plouf.samples[(aa+2):(aa+1+nimp)]
                                   imps <- as.integer(unlist(strsplit(imp,split=':')))
                                   ns <- imps[seq(2,2*nimp,2)]
                                   ds <- imps[seq(1,2*nimp,2)]
                                }
                                return(list(ns=ns, ds=ds)) })
      }

      # a) left and right reads overlap or stop/start into same exon or neighbour exons -- EASY case
      if(lastleft>=firstright-1){
         vecbin <- unique(sort(c(oneleft,oneright)))
         newbin <- paste(vecbin, collapse='-')
         longtype[[newbin]] <- vecbin
         numbertype[[newbin]] <- c(numbertype[[newbin]],npair)
         if(!is.null(samples.name)){
            if(is.null(numbertype.samples[[newbin]])) { numbertype.samples[[newbin]] <- rep(0, n.samples) }
            numbertype.samples[[newbin]] <- numbertype.samples[[newbin]] + sapply(info.samples, FUN=function(ll) sum(ll$ns)) # count of the bin for each sample
         }
      }
      
      # b) left and right reads stop/start not into neighbours bins -- MIGHT BE AMBIGUOUS
      if(firstright>lastleft+1){
         # b1) only jump an intron (active when annotation in input)
         if(sum(count.exons[(lastleft+1):(firstright-1)])==0){  # be less stringent ????
            vecbin <- unique(sort(c(oneleft,oneright)))
            newbin <- paste(vecbin, collapse='-')
            longtype[[newbin]] <- vecbin
            numbertype[[newbin]] <- c(numbertype[[newbin]],npair)
            if(!is.null(samples.name)){
               if(is.null(numbertype.samples[[newbin]])) { numbertype.samples[[newbin]] <- rep(0, n.samples) }
               numbertype.samples[[newbin]] <- numbertype.samples[[newbin]] + sapply(info.samples, FUN=function(ll) sum(ll$ns))
            }
         }
         # b2) complicated cases
         if(sum(count.exons[(lastleft+1):(firstright-1)])>0){
            djump <- alldist - (tophat.exons[firstright,1] - tophat.exons[lastleft,2]) # effective insert when jump all inside bins
            ppp <- which(count.exons[(lastleft+1):(firstright-1)]>0)
            addno <- sum(len.exons[(lastleft+1):(firstright-1)][ppp])
            dnojump <- djump + addno # effective insert when no jump at all (exept introns of course)   # be less stringent ???
            seuil.nojump <- min(gap-std, gap+std-addno)
            ind.nojump <- (djump<seuil.nojump) # indexes of no jumping reads

            if(any(ind.nojump)){
               vecbin <- unique(sort(c(oneleft,((lastleft+1):(firstright-1))[ppp],oneright)))
               newbin <- paste(vecbin, collapse='-')
               longtype[[newbin]] <- vecbin
               numbertype[[newbin]] <- c(numbertype[[newbin]],npair[ind.nojump])
               if(!is.null(samples.name)){
               sumsamp <- sapply(info.samples, FUN=function(ll){
                      rsum <- 0
                      if(!is.null(ll$ds)) {
                         tt <- sapply(ll$ds, FUN=function(dd) ind.nojump[alldist==dd])
                         if(any(tt)) rsum <- sum(ll$ns[tt])
                      }
                      return(rsum) })
                  if(is.null(numbertype.samples[[newbin]])) { numbertype.samples[[newbin]] <- rep(0, n.samples) }
                  numbertype.samples[[newbin]] <- numbertype.samples[[newbin]] + sumsamp
               }
            }
            
            seuil.jump <- max(gap-std, gap+std-addno)
            ind.jump <- (djump>seuil.jump) # indexes of jumping reads
            if(any(ind.jump)){
               vecbin <- unique(sort(c(oneleft,oneright)))
               newbin <- paste(vecbin, collapse='-')
               longtype[[newbin]] <- vecbin
               numbertype[[newbin]] <- c(numbertype[[newbin]],npair[ind.jump])
               if(!is.null(samples.name)){
               sumsamp <- sapply(info.samples, FUN=function(ll){
                      rsum <- 0
                      if(!is.null(ll$ds)) {
                         tt <- sapply(ll$ds, FUN=function(dd) ind.jump[alldist==dd])
                         if(any(tt)) rsum <- sum(ll$ns[tt])
                      }
                      return(rsum) })
                  if(is.null(numbertype.samples[[newbin]])) { numbertype.samples[[newbin]] <- rep(0, n.samples) }
                  numbertype.samples[[newbin]] <- numbertype.samples[[newbin]] + sumsamp
               }
            }
         }
      }
   }

   # II) second pass - readjust bin
   # if a type if strictly inside another type then transfer it
   if(length(longtype)>0){

      longlength <- sapply(longtype, length)
      inds <- sort(longlength, index.r=T)$ix
      numbersum <- sapply(numbertype, sum)

      toremove <- rep(FALSE, length(longtype))
      toadd <- integer(length(numbertype))

      if(!is.null(samples.name)){
         numbersum.samples <- lapply( 1:n.samples, FUN=function(ii) sapply(numbertype.samples, FUN=function(ss) ss[ii], USE.NAMES=F))
         toadd.samples <- lapply(samples.name, FUN=function(ll) return(toadd))
      }

      for(jj in 1:length(inds)){
         ii <- inds[jj]
         inbin <- longtype[[ii]]
         if(inbin[1]>1 & inbin[length(inbin)]<n.exons){ # cannot be inside if start by 1st exon or and by last exon
            totest <- inds[jj:length(inds)]
            # select bins such that inbin is strictly inside
            test.inside <- sapply(longtype[totest], FUN=function(outbin) bin.inside(inside.bin=inbin, outside.bin=outbin))
            # check
            if(any(test.inside)){
               toremove[ii] <- TRUE # the bin will be transfered
               totransfer <- numbersum[ii]
               testpos0 <- totest[which(test.inside)]
               testpos <- testpos0[which.min( longlength[testpos0] )]
               proportion <- numbersum[testpos]/sum(numbersum[testpos])
               toadd[testpos] <- proportion*totransfer
               if(!is.null(samples.name)){
                  # do it for each sample
                  for(nn in 1:n.samples){
                     tot <- numbersum.samples[[nn]][ii]
                     norms <- sum(numbersum.samples[[nn]][testpos])
                     if(norms==0) {
                        prop <- proportion # if no read at all in the sample, use global proportion
                     } else { 
                        prop <- numbersum.samples[[nn]][testpos]/norms # else use sample proportions
                     }
                     toadd.samples[[nn]][testpos] <- prop*tot
                  }
               }
            }
         }
      }

      keeplolo <- longtype[!toremove]
      count.first <- (numbersum + toadd)[!toremove]
      if(!is.null(samples.name)){
         count.first.samples <- lapply(1:n.samples, FUN=function(nn) (numbersum.samples[[nn]] + toadd.samples[[nn]])[!toremove])
      }

   }

   return(list(longtype=longtype, keeplolo=keeplolo, count.first=count.first, count.first.samples=count.first.samples))

}
