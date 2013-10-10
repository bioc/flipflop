##========== Splicing graph for single-end reads ===========##

#*** Create all the 'bins' i.e all the node of the splicing graph
#    and the junctions between bins ***
# INPUT:   - observed bins
# OUTPUT:  - all bins (oberved + created)
#          - splicing graph
#          - vectors of count and cost on the bins

graph_single <- function(binnames, goodtype, readlen, len.exons, n.exons, tophat.exons, use_TSSPAS){
   binary <- goodtype[,1:n.exons, drop=F]
   binreste <- binnames # observed bins in the initial queue
   allbins <- binnames # all bins - observed and created - may be completed
   list.junctions <- vector() # directed junctions between bins

   while(length(binreste)>0){ # Look for appropriate junction between bins
      currentbin <- binreste[1]
      attrape <- as.numeric(strsplit(currentbin, split='-')[[1]])
      debut <- attrape[1]
      fin <- attrape[length(attrape)]
      if(length(attrape)>1){ # multi-exons type
         rightbin <- paste(attrape[-1], collapse='-') # Look for shorter going out junctions
         lright <- sum(len.exons[attrape[-1]])
         if(lright >= readlen){ # rightbin is compatible
            if(match(rightbin, allbins, nomatch=0)==0){ # rightbin not in the observed set
               allbins <- c(allbins, rightbin)  # add the shorter bin
               binreste <- c(binreste, rightbin) 
            } 
            list.junctions <- rbind(list.junctions, c(currentbin, rightbin)) # add the junction
         }
         if(lright<readlen & fin<n.exons){ # rightbin is not compatible
            nextright <- vector()
            for(qq in 1:length(attrape)){ # test all possible suffix
               suffix <- attrape[qq:length(attrape)]
               lala <- apply(matrix(binary, ncol=n.exons), 1, FUN=function(xx)  find.prefix(suffix,xx, tophat.exons)) # isolate same prefix
               if(sum(lala)>0){
                  nextright <- c(nextright,apply(matrix(binary[lala,], ncol=n.exons), 1, FUN=function(yy) min(setdiff(which(yy==1),suffix)))) # first non zero bin after suffix
               }
            }
            nextright <- unique(nextright)
            if(length(nextright)>0){
               if(lright==(readlen-1)){ # particular case
                  for(addbin in nextright){
                     longerbin <- paste(c(rightbin,addbin),collapse='-')
                     if(match(longerbin, allbins,nomatch=0)==0){ # add the longer bin
                        allbins <- c(allbins, longerbin)
                        binreste <- c(binreste, longerbin)
                     }
                     list.junctions <- rbind(list.junctions, c(currentbin, longerbin)) # add the junction
                  }
               }
               if(lright<(readlen-1)){ # Look for longer going out junction
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
      binreste <- binreste[binreste!=currentbin] # remove currentbin to queue
   }
   # Starting bins with no entering junctions - look if there is a long enough prefix to add a coming in junction
   if(length(list.junctions)>0){
      noentrant <- allbins[match(allbins, list.junctions[,2], nomatch=0)==0]
   }
   if(length(list.junctions)==0){
      noentrant <- allbins
   }
   for(firstbin in noentrant){
      attrape <- as.numeric(strsplit(firstbin, split='-')[[1]])
      lleft <- sum(len.exons[attrape[-length(attrape)]])
      gobin <- firstbin
      while(length(attrape)>1 & lleft >=readlen){
         leftbin <- paste(attrape[-length(attrape)], collapse='-')
         allbins <- c(allbins, leftbin) # add leftbins
         list.junctions <- rbind(list.junctions, c(leftbin, gobin)) # add the junction
         attrape <- attrape[-length(attrape)]
         lleft <- sum(len.exons[attrape[-length(attrape)]])
         gobin <- leftbin
      }
   }

   ## Adjacency matrix
   n.nodes <- length(allbins)
   CC <- diag(0, n.nodes)
   rownames(CC) <- colnames(CC) <- allbins
   CC[list.junctions] <- 1
   CC <- t(CC)
   sp.Adj <- as(CC,'CsparseMatrix')
   ## Counts on bins
   count <- matrix(0,nrow=n.nodes,ncol=1)
   rownames(count) <- allbins
   count[binnames,] <- goodtype[binnames,n.exons+1]
   ## Cost on bins
   len.bin <- sapply(allbins, FUN=function(zz) calculate.cost(len.exons, zz, readlen))
   len <- matrix(len.bin, nrow=n.nodes, ncol=1)
   rownames(len) <- allbins
   ## Weights for START and STOP:
   wg.start <- wg.stop <- rep(Inf, n.nodes)
   toto <- sapply(allbins, FUN=function(tt) startstop(len.exons,tt,readlen,count))
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


