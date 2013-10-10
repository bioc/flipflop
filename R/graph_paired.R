##========== Splicing graph for paired-end reads ===========##

#*** Create all the 'bins' i.e all the node of the splicing graph
#    and the junctions between bins ***
# INPUT:   - observed bins
# OUTPUT:  - all bins (oberved + created)
#          - splicing graph
#          - vectors of count and cost on the bins


graph_paired <- function(binnames, len.exons, n.exons, tophat.exons, binary, frag, numbertype, use_TSSPAS){
   list.junctions <- vector()
   binreste <- binnames # queue
   allbins <- binnames # all bins - observed and created -

   while(length(binreste)>0){ # Look for appropriate junction between nodes
      currentbin <- binreste[1]
      attrape <- as.numeric(strsplit(currentbin, split='-')[[1]])
      debut <- attrape[1]
      fin <- attrape[length(attrape)]
      if(length(attrape)>1){ # multi-exons type
         rightbin <- paste(attrape[-1], collapse='-') # Look for shorter going out junctions
         lright <- sum(len.exons[attrape[-1]])
         if(match(rightbin, allbins, nomatch=0)>0){ # rightbin is in the observed set    
            list.junctions <- rbind(list.junctions, c(currentbin, rightbin)) # add the junction
         }
         if(match(rightbin, allbins, nomatch=0)==0 & lright>frag){ # rightbin not in the observed set and considered as compatible
            allbins <- c(allbins, rightbin) # add the bin
            binreste <- c(binreste, rightbin)
            list.junctions <- rbind(list.junctions, c(currentbin, rightbin))
         }
         if(fin<n.exons){ 
            nextright <- vector()
            for(qq in 1:length(attrape)){ # test all possible suffix of currentbin
               suffix <- attrape[qq:length(attrape)]
               lala <- apply(matrix(binary, ncol=n.exons), 1, FUN=function(xx)  find.prefix(suffix,xx, tophat.exons)) # isolate same prefix
               if(sum(lala)>0){
                  nextright <- c(nextright,apply(matrix(binary[lala,], ncol=n.exons), 1, FUN=function(yy) min(setdiff(which(yy==1),suffix)))) # first non zero bin after suffix
               }
            }
            nextright <- unique(nextright)
            if(length(nextright)>0){
               for(addbin in nextright){
                  longerbin <- paste(c(currentbin,addbin),collapse='-')
                  if(match(longerbin, allbins,nomatch=0)>0){ # longer bin is in the observed set
                     list.junctions <- rbind(list.junctions, c(currentbin, longerbin))
                  }
                  if(match(longerbin, allbins,nomatch=0)==0 & match(rightbin, allbins,nomatch=0)==0){ # longer bin is not in the observed set and neither rightbin
                     allbins <- c(allbins, longerbin)
                     binreste <- c(binreste, longerbin)
                     list.junctions <- rbind(list.junctions, c(currentbin, longerbin))
                  }
               }
            }
         }
      }
      if(length(attrape)==1){ # mono-exon type
         if(debut<n.exons){ # Look for longer going out junction
            lala <- apply(matrix(binary, ncol=n.exons), 1, FUN=function(xx)  find.prefix(debut,xx,tophat.exons)) # isolate same prefix
            if(sum(lala)>0){
               nextright <- apply(matrix(binary[lala,], ncol=n.exons), 1, FUN=function(yy) min(setdiff(which(yy==1),debut))) # first non zero bin after suffix
               nextright <- unique(nextright)
               for(addbin in nextright){
                  longerbin <- paste(c(currentbin,addbin),collapse='-')
                  if(match(longerbin, allbins,nomatch=0)==0){ # add the longer bin
                     allbins <- c(allbins, longerbin)
                     binreste <- c(binreste, longerbin)
                  }
                  list.junctions <- rbind(list.junctions, c(currentbin, longerbin)) # add the junction
               }
            }
         }
      }
      binreste <- binreste[binreste!=currentbin] # remove currentbin to queue  (-1 suffirait non ?)
   }
   ## PAIRED-END : je dois faire un prefix search pour trouver les éventuels jonctions entrantes.
   # Starting bins with no entering junctions - look if there is a long enough prefix to add a coming in junction
   if(length(list.junctions)>0){
      noentrant <- allbins[match(allbins, list.junctions[,2], nomatch=0)==0]
   }
   if(length(list.junctions)==0){
      noentrant <- allbins
   }
   while(length(noentrant)>0){
      firstbin <- noentrant[1]
      attrape <- as.numeric(strsplit(firstbin, split='-')[[1]])
      if(length(attrape)>1){
         allprefix <- vector(mode='character',length=(length(attrape)-1))
         for(rr in 1:(length(attrape)-1)){
            allprefix[rr] <- paste(attrape[1:rr], collapse='-')
         }
         prefix.match <- match(allprefix, allbins, nomatch=0)
         if(sum(prefix.match)>0){ # one of prefix is in the existing type set
            newentrant <- allprefix[length(allprefix)]
            list.junctions <- rbind(list.junctions, c(newentrant, firstbin)) # add the junction 
            if(prefix.match[length(prefix.match)]==0){ # the prefix for new junction do not exist 
               allbins <- c(allbins, newentrant)
               noentrant <- c(noentrant, newentrant)
            }
         }
         if(sum(prefix.match)==0 & sum(len.exons[attrape[-length(attrape)]])>frag){ # none of the prefix exists but the longest one is compatible
            newentrant <- allprefix[length(allprefix)]
            list.junctions <- rbind(list.junctions, c(newentrant, firstbin))
            allbins <- c(allbins, newentrant)
            noentrant <- c(noentrant, newentrant)
         }
      }
      noentrant <- noentrant[-1] # pop
   }
   list.junctions <- unique(list.junctions) # peut faire mieux ?

   ## Adjacency matrix
   n.nodes <- length(allbins)
   CC <- diag(0, n.nodes)
   rownames(CC) <- colnames(CC) <- allbins
   CC[list.junctions] <- 1
   CC <- t(CC)
   sp.Adj <- as(CC,'CsparseMatrix')
   ## Count on bins
   count <- matrix(0,nrow=n.nodes,ncol=1)
   rownames(count) <- allbins
   count[binnames,] <- sapply(numbertype,sum)
   ## Cost on bins
   len.init <- sapply(allbins, FUN=function(zz) calculate.cost(len.exons, zz, frag))
   len <- matrix(len.init, nrow=n.nodes, ncol=1)
   rownames(len) <- allbins
   ## Weights for START and STOP:
   wg.start <- wg.stop <- rep(Inf, n.nodes)
   toto <- sapply(allbins, FUN=function(tt) startstop_paired(tt,allbins)) 
   allstart <- allbins[toto[1,]]
   allstop <- allbins[toto[2,]]
   if(use_TSSPAS==0){
      wg.start[allbins %in% allstart] <- 1
      wg.stop[allbins %in% allstop] <- 0
   }
   if(use_TSSPAS==1){
      ## TSS pas de jonction entrante (1 ligne de 0)
      ind.sometss <- apply(CC,1,sum)==0 
      sometss <- allbins[ind.sometss]
      ## PAS pas de jonction sortante (1 col de 0)
      ind.somepas <- apply(CC,2,sum)==0
      somepas <- allbins[ind.somepas]
      premier.tss <- sapply(sometss, FUN=function(v) firstexon(v))
      dernier.pas <- sapply(somepas, FUN=function(v) lastexon(v)) 
      alltss <- allstart[sapply(allstart, FUN=function(u) findtss(u,premier.tss, tophat.exons, len.exons))]
      allpas <- allstop[sapply(allstop, FUN=function(u) findpas(u,dernier.pas, tophat.exons, len.exons))]
      wg.start[allbins %in% alltss] <- 1
      wg.stop[allbins %in% allpas] <- 0
   }
   graph <- list()
   graph[['weights']] <- sp.Adj*1e-24
   graph[['start_weights']] <- wg.start
   graph[['stop_weights']] <- wg.stop

   return(list(allbins=allbins, count=count, len=len, graph=graph))
}
