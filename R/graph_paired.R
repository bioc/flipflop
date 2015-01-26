# --------------- Splicing graph for paired-end reads --------------- #

# *** Create all the 'bins' i.e all the nodes of the splicing graph and the junctions between bins ***
#
# INPUT:   - observed bins: binlist
#          - count on the observed bins: count.first
#          - exons: length len.exons and positions tophat.exons
#          - mean fragment length: frag
#          - TSS and PAS from the annotations
# OUTPUT:  - allbins (oberved + created)
#          - splicing graph
#          - vectors of count and cost on the bins: count and len


# ( 13/01/14 new version: manipulate clever structure
# do not use multiple rbind anymore
# do not use multiple strsplit / paste
# change paradigm: map with an entry for each exon number giving the list of type that starts by that exon )


graph_paired <- function(binlist, count.first, frag, len.exons, n.exons, tophat.exons, use_TSSPAS, TSSref, PASref) {

   # Look for appropriate junctions between bins
   # And potentially create new bins

   # create a list of list
   first.exons <- sapply(binlist, FUN=function(bb) bb[1])
   unique.first.exons <- unique(first.exons)
   list.byexons <- vector(mode='list', length=n.exons)
   list.junctions <- list()
   list.byexons[unique.first.exons] <- lapply( unique.first.exons, FUN=function(ee) binlist[first.exons==ee] )

   length.all <- sapply(list.byexons, length)
   length.done <- rep(0, n.exons)
   length.init <- length.all # !! 

   list.byexons00 <- list.byexons

   njunc <- 0

   while(sum(length.all - length.done)>0) {
  
      debut <- which( length.all - length.done !=0 )[1]
      indcurrent <- length.done[debut]+1
      attrape <- list.byexons[[debut]][[indcurrent]]
      fin <- attrape[length(attrape)]

      #-----------------
      # multi-exons bin:
      #-----------------

      if(length(attrape)>1){

         rightbin <- attrape[-1] # Look for shorter going out junctions
         lright <- sum(len.exons[rightbin])

         rightbin.first <- rightbin[1]
         testmatch <- quickmatch(rightbin,list.byexons[[rightbin.first]])

         # rightbin is in the observed set:
         if(testmatch){
            njunc <- njunc+1 
            list.junctions[[njunc]] <- c(debut,indcurrent,rightbin.first,testmatch)
         }

         # rightbin not in the observed set and considered as compatible:
         if(!testmatch & lright>frag){
            length.all[rightbin.first] <- length.all[rightbin.first]+1
            testmatch <- length.all[rightbin.first]
            # RESOUDRE LE PROBLEME DEBILE DE R QUI NE FAIT PAS UNE LISTE QUAND 1er ELEMENT de taille 1 !!
            list.byexons[[rightbin.first]][[testmatch]] <- integer(2)
            list.byexons[[rightbin.first]][[testmatch]] <- rightbin
            njunc <- njunc+1 
            list.junctions[[njunc]] <- c(debut,indcurrent,rightbin.first,testmatch)
         }
         testmatch.right <- testmatch
         
         # suffix/prefix search:
         if(fin<n.exons){
            nextright <- vector()
            for(qq in 1:length(attrape)){ # test all possible suffix of currentbin
               suffix <- attrape[qq:length(attrape)]
               if(length(list.byexons00[[suffix[1]]])>0){ # consider non empty key 
                  lala <- sapply(list.byexons00[[suffix[1]]], FUN=function(bin)  find.prefix(suffix,bin,tophat.exons)) # isolate same prefix
                  if(sum(lala)>0){
                     nextright <- c(nextright,lala[lala>0]) # first non zero bin after suffix
                  }
               }
            }
            # check if the exon right after is in the list
            # if not AND consecutive genomic postion then add it # TODO seulement no annot
            if( match(fin+1, nextright,nomatch=0)==0 & 
               length(list.byexons00[[fin+1]])>0 &
               tophat.exons[fin,2]==tophat.exons[fin+1,1] ){
               nextright <- c(nextright, fin+1)
            }
            nextright <- as.integer(unique(nextright))
            all.longerbin <- lapply(nextright, FUN=function(nn) c(attrape,nn))
            for(longerbin in all.longerbin){
               longfirst <- longerbin[1]
               testmatch <- quickmatch(longerbin,list.byexons[[longfirst]])
               if(testmatch){ # longerbin is in the observed set
                  njunc <- njunc+1
                  list.junctions[[njunc]] <- c(debut,indcurrent,longfirst,testmatch)
               }
               if(!testmatch & !testmatch.right){ # longerbin is not in the observed set AND neither rightbin
                  length.all[longfirst] <- length.all[longfirst]+1
                  testmatch <- length.all[longfirst]
                  list.byexons[[longfirst]][[testmatch]] <- longerbin
                  njunc <- njunc+1
                  list.junctions[[njunc]] <- c(debut,indcurrent,longfirst,testmatch)
               }
            }
         }
      }

      #-----------------
      # mono-exon bin:
      #-----------------

      if(length(attrape)==1){ # mono-exon type
         if(debut<n.exons){ # Look for longer going out junction
            nextright <- vector()
            if(length(list.byexons00[[debut]])>0){
               lala <- sapply(list.byexons00[[debut]], FUN=function(bin)  find.prefix(debut,bin,tophat.exons)) # isolate same prefix
               if(sum(lala)>0){
                  nextright <- lala[lala>0] # first non zero bin after suffix
               }
            }
            # check if the exon right after is in the list
            # if not AND consecutive genomic postion then add it
            if(match(fin+1, nextright,nomatch=0)==0 & 
               length(list.byexons00[[fin+1]])>0 &
               tophat.exons[fin,2]==tophat.exons[fin+1,1] ){
               nextright <- c(nextright, fin+1)
            }
            nextright <- as.integer(unique(nextright))
            all.longerbin <- lapply(nextright, FUN=function(nn) c(attrape,nn))
            for(longerbin in all.longerbin){
               testmatch <- quickmatch(longerbin,list.byexons[[debut]])
               if(!testmatch){ # add the longer bin
                  length.all[debut] <- length.all[debut]+1
                  testmatch <- length.all[debut]
                  list.byexons[[debut]][[testmatch]] <- longerbin
               }
               njunc <- njunc+1
               list.junctions[[njunc]] <- c(debut,indcurrent,debut,testmatch)
            }
         }
      }
      length.done[debut] <- length.done[debut]+1
   }
   

   ## PAIRED-END : je dois faire un prefix search pour trouver les éventuels jonctions entrantes.
   # Starting bins with no entering junctions - look if there is a long enough prefix to add a coming in junction
   for(debut in 1:n.exons){
      res <- sapply(list.junctions, FUN=function(lj) {
                    if(lj[3]==debut){
                       return(lj[4])
                    } else {
                       return(0)
               } } )
      noentrant <- setdiff(1:length(list.byexons[[debut]]), res)
      for(ind in noentrant){
         ind0 <- ind
         attrape <- list.byexons[[debut]][[ind]]
         if(length(attrape)>1){
            gobin <- attrape[-length(attrape)]
            lleft <- sum(len.exons[gobin])
            # consider all prefix of "attrape"
            allprefix <- vector(mode='list', length=(length(attrape)-1))
            for(rr in 1:(length(attrape)-1)){
               allprefix[[rr]] <- attrape[1:rr]
            }
            alltest <- sapply(allprefix, FUN=function(ppb) quickmatch(ppb,list.byexons[[debut]]) )
            # one of the prefix is in the existing set OR the longest one is compatible:
            while( lleft >=frag | sum(alltest)>0 ){
               testmatch <- alltest[length(alltest)]
               if(!testmatch){ # the prefix does not already exist -- add bin
                  length.all[debut] <- length.all[debut]+1
                  testmatch <- length.all[debut]
                  list.byexons[[debut]][[testmatch]] <- gobin
               }
               njunc <- njunc+1
               list.junctions[[njunc]] <- c(debut,testmatch,debut,ind0)
               ind0 <- testmatch
               attrape <- attrape[-length(attrape)]
               gobin <- attrape[-length(attrape)]
               lleft <- sum(len.exons[gobin])
               alltest <- alltest[-length(alltest)]
            }
         }
      }
   }
   
   ## All bins in a list
   all.list <- unlist(list.byexons, recursive=F, use.names=F)

   ## Adjacency matrix  # faire directement une matrice sparse !!
   n.nodes <- length(all.list)
   if(njunc>0){
      pos.junc <- sapply(list.junctions, FUN=function(lj){
                         if(lj[1]>1){
                            pos.first <- sum( length.all[1:(lj[1]-1)] ) + lj[2]
                         } else {
                            pos.first <- lj[2]
                         }
                         if(lj[3]>1){
                            pos.second <- sum( length.all[1:(lj[3]-1)] ) + lj[4]
                         } else {
                            pos.second <- lj[4]
                         }
                         return(c(pos.second, pos.first))
               })

      sp.Adj <- sparseMatrix(i=pos.junc[1,],j=pos.junc[2,], dims=c(n.nodes,n.nodes), giveCsparse=T)
   } else {
      sp.Adj <- sparseMatrix(i=NULL,j=NULL, dims=c(n.nodes,n.nodes), giveCsparse=T)
   }

   allbins <- sapply(all.list, FUN=function(aa) paste(aa, collapse='-'))
   rownames(sp.Adj) <- allbins

   ## Count on bins  # le faire directement dans la boucle ?
   count <- matrix(0,nrow=n.nodes,ncol=1)
   indcount <-  sapply(1:length(first.exons), FUN=function(ii){
                       debut <- first.exons[ii]
                       pos <- sum(first.exons[1:ii]==debut)
                       if(debut==1){
                          return( pos )
                       } else {
                          return( sum(length.all[1:(debut-1)]) + pos )
                       }
               })

   count[indcount] <- count.first
   rownames(count) <- allbins

   ## Cost on bins
   len.bin <- sapply(all.list, FUN=function(zz) calculate.cost(len.exons, zz, min(frag,sum(len.exons)))) # 14/01/15 not only zero weights
   len <- matrix(len.bin, nrow=n.nodes, ncol=1)
   rownames(len) <- allbins 
   
   ## Weights for START and STOP:
   wg.start <- wg.stop <- rep(Inf, n.nodes)
   toto <- sapply(all.list, FUN=function(tt) startstop_paired(tt,list.byexons)) 
   indalltss <- which(toto[1,]) # index of possible TSS bins
   indallpas <- which(toto[2,]) # index of possible PAS bins
   if(use_TSSPAS==0){
      wg.start[indalltss] <- 1
      wg.stop[indallpas] <- 0
   }
   if(use_TSSPAS==1){
      ## TSS no entering junctions (1 line of 0)
      ind.sometss <- which(!rowSums(sp.Adj))
      ## PAS no outgoing junctions (1 col of 0)
      ind.somepas <- which(!colSums(sp.Adj))
      if(is.null(TSSref) | is.null(PASref)) { # no annotation
         # add direct neightbours of previously selected TSS PAS 
         exon.tss <- sapply(ind.sometss, FUN=function(ii) firstexon(all.list[[ii]]))
         exon.pas <- sapply(ind.somepas, FUN=function(ii) lastexon(all.list[[ii]])) 
         indtss <- indalltss[sapply(indalltss, FUN=function(ii) findtss(all.list[[ii]],exon.tss, tophat.exons, len.exons))]
         indpas <- indallpas[sapply(indallpas, FUN=function(ii) findpas(all.list[[ii]],exon.pas, tophat.exons, len.exons))]
      } else { # use known TSS PAS from annotation
         # add the known TSS PAS
         indtss <- union(ind.sometss,indalltss[sapply(indalltss, FUN=function(ii) all.list[[ii]][1] %in% TSSref)])
         indpas <- union(ind.somepas, indallpas[sapply(indallpas, FUN=function(ii) tail(all.list[[ii]],1) %in% PASref)])
      }
      wg.start[indtss] <- 1
      wg.stop[indpas] <- 0
   }

   graph <- vector(mode="list", length=3)
   names(graph) <- c('weights','start_weights','stop_weights')
   graph[['weights']] <- sp.Adj*1e-24
   graph[['start_weights']] <- wg.start
   graph[['stop_weights']] <- wg.stop

   return(list(allbins=allbins, count=count, len=len, graph=graph))
}
