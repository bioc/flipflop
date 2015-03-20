
WrapperFlipFlop <- function(
                     inpf, # the input open file
                     outf, # the output open file
                     outf.fpkm,
                     outf.count,
                     output.type,
                     NN,
                     samples.name,
                     n.samples,
                     ind.samples,
                     mergerefit,
                     paired,
                     frag,
                     std,                                                      
                     minReadNum,
                     minFragNum,
                     mc.cores,
                     verbose,
                     verbosepath,
                     cutoff,
                     BICcst,
                     delta,
                     use_TSSPAS,
                     max_isoforms
                     ) {

      ##=======================================
      ## Set up parameters
      ##=======================================
      #delta <- 1e-5 # Smoothing parameter for the Poisson Model # now an option
      param <- list()
      param$mode_decomposition <- 3 # Flow decomposition mode (1, 2 or 3)
      param$regul <- 'graph-path-conv2'
      param$verbose <- FALSE
      param$loss <- 'poisson-weighted'
      param$delta <- delta
      param$pos <- TRUE
      param$tol <- 1e-10
      ##=======================================

      beta.fpkm <- list()
      beta.expcount <- list()
      timer.all <- vector()
      granges <- GenomicRangesList()
      mm <- 0
      scan(inpf, what=character(0), nlines=3, quiet=TRUE)

      while(1){
         mm <- mm+1

         ## Name-Instance
         name <- scan(inpf, what=character(0), nlines=1, quiet=TRUE)
         print(name)
         nom <- paste('Inst', name[2], sep='')

         if(length(name)==0) break

         ## Boundary
         bound <- scan(inpf, what=character(0), nlines=1, quiet=TRUE)
         chr <- bound[2]
         strand <- bound[5]

         ## Read Lenght
         readinfo <- scan(inpf, what=character(0), nlines=1, quiet=TRUE)
         readlen <- as.integer(readinfo[2])
         gap <- frag - 2*readlen

         ## Number of (sub)-exons
         segs <- scan(inpf, what=character(0), nlines=1, quiet=TRUE)
         n.exons <- as.integer(segs[2])
         print(paste('sub-exons', n.exons, sep=':'))

         ## Exons Position
         exonsinfo <- scan(inpf, what=numeric(0), nlines=n.exons, quiet=TRUE)
         dim(exonsinfo) <- c(9,n.exons)
         exonsinfo <- t(exonsinfo)
         tophat.exons <- exonsinfo[1:n.exons,1:2, drop=FALSE]
         tophat.exons[,2] <- tophat.exons[,2] + 1

         ## Count and Length on Exons
         len.exons <- exonsinfo[1:n.exons,3]
         count.exons <- exonsinfo[1:n.exons,4]

         ## Annotated References
         ref <- scan(inpf, what=character(0), nlines=1, quiet=TRUE)
         nref <- as.integer(ref[2])
         TSSref <- PASref <- NULL
         if(nref>=1){
            # so far we don't use reference information for more than exon boundaries and (optionnaly) TSS and PAS
            listref <- scan(inpf, what=character(0), nlines=nref, quiet=TRUE) 
            dim(listref) <- c((n.exons+2),nref)
            listref <- listref[1:n.exons,,drop=F]
            TSSref <- unique(apply(listref, 2, FUN=function(v) which(v==1)[1]))
            PASref <- unique(apply(listref, 2, FUN=function(v) tail(which(v==1),1)))
         }

         ## Total number of reads for this gene
         nreads <- scan(inpf, what=character(0), nlines=1, quiet=TRUE)
         nreads <- as.integer(nreads[2])

         ## Read type
         type <- scan(inpf, what=character(0), nlines=1, quiet=TRUE)
         ntype <- as.integer(type[2])
         if(ntype > 0){
            readtype_all <- scan(inpf, what=character(0), nlines=ntype, quiet=TRUE)
            dim(readtype_all) <- c(length(readtype_all)/ntype, ntype)
            readtype_all <- t(readtype_all)
            readtype <- readtype_all[,1:(n.exons+1), drop=F]
            readtype <- apply(readtype, 2, as.double)  # readtype contains all types of reads for the sum of samples (as.double otherwise bug)
            if(!is.matrix(readtype)) readtype <- matrix(readtype, ncol=length(readtype))
            n.samples.max <- ncol(readtype_all)-ncol(readtype)-1 # the total number of samples in the provied instance file (potentially more than the given samples)
         } else {
            VALID <- FALSE
         }

         ## Paired-End type
         petype <- scan(inpf, what=character(0), nlines=1, quiet=TRUE)
         npereads <- as.integer(petype[2])
         npetype <- as.integer(petype[3])
         # Skip instance when not enough mapped reads or mapped fragments 
         # (but don't skip if your reads are supposed to be paired but there is more than 1000 single-end reads)
         if(nreads<minReadNum | (paired==TRUE & npereads<minFragNum & nreads<1000)){
            if(npetype>0){
               scan(inpf, what=character(0), nlines=((2+n.samples.max)*npetype), quiet=TRUE) # Skip paired-end type information
            }
            next
         }
         # Consider reads as single-end (npetype==0 means there is no pair and for npetype>300 the paired graph creation is too long):
         eff.paired <- FALSE
         if(npetype>0){
            if(paired==FALSE | npetype>=300) {
               scan(inpf, what=character(0), nlines=((2+n.samples.max)*npetype), quiet=TRUE) # Skip paired-end type information
            } else {
               # Paired-end reads:
               eff.paired <- TRUE
               # Construct the bins (consider pairs as 'long single' with some heuristics if necessary)
               bins.res <- bins_paired(inpf=inpf, npetype=npetype, tophat.exons=tophat.exons, count.exons=count.exons, 
                                       len.exons=len.exons, n.exons=n.exons, readtype=readtype, gap=gap, std=std, 
                                       n.samples.max=n.samples.max, samples.name=samples.name, n.samples=n.samples)
               longtype <- bins.res$longtype
               VALID <- (length(longtype)>0)
               keeplolo <- bins.res$keeplolo
               count.first <- bins.res$count.first
            }
         }

         tic()
         ##==========================================
         ## Create Splicing Graph
         ## NODE:=READTYPE:=BIN
         ##==========================================
         if(!eff.paired & nreads>0){
            lolo <- lapply(1:nrow(readtype), FUN=function(v) which(readtype[v,1:n.exons]==1))
            ## sanity check 2 --- valid bins ----------------------------------
            validbins <- sapply(lolo, FUN=function(v){
                                if(length(v)<=2){ return((sum(len.exons[v])>=readlen)) }
                                if(length(v)>=3){ vv <- v[2:(length(v)-1)] ; return(((sum(len.exons[v])>=readlen) & (sum(len.exons[vv])<(readlen-1)))) }
                                       })
            VALID <- (sum(validbins)>0)
         }

         if(!VALID) { # don't go if you do not have reads or any bins
            path.select <- matrix(0, nrow=n.exons, ncol=1)
            beta.select <- beta.raw <- rep(0, n.samples)
         } else {
            # Create bins and graph
            if(!eff.paired) {
               keeplolo <- lolo[validbins]
               count.first <- readtype[validbins,n.exons+1]
               graph_creation <- graph_single(binlist=keeplolo, count.first=count.first, readlen=readlen,
                                              len.exons=len.exons, n.exons=n.exons, tophat.exons=tophat.exons,
                                              use_TSSPAS=use_TSSPAS, TSSref=TSSref, PASref=PASref)
            } else {
               graph_creation <- graph_paired(binlist=keeplolo, count.first=count.first, frag=frag, 
                                              len.exons=len.exons, n.exons=n.exons, tophat.exons=tophat.exons,
                                              use_TSSPAS=use_TSSPAS, TSSref=TSSref, PASref=PASref)
            }
            allbins <- graph_creation$allbins
            n.nodes=length(allbins)
            count <- graph_creation$count
            len <- graph_creation$len
            graph <- graph_creation$graph
            
            # save matrix of count for each samples ---------------------
            if(!is.null(samples.name)){
               indcount <- graph_creation$indcount
               count.samples <- matrix(0, nrow=n.nodes, ncol=n.samples)
               if(!eff.paired) {
                  goodtype_all <- readtype_all[validbins,(n.exons+2):(n.exons+1+n.samples.max),drop=F]
                  count.first.samples <- apply(goodtype_all, 1, FUN=function(v){
                                               pp <- strsplit(v, split='=')
                                               tt <- sapply(pp, FUN=function(k) as.double(k[2]))  # double ?? 
                                               return(tt)})
                  count.samples[indcount,] <- t(count.first.samples)[,ind.samples]
               } else {
                  count.first.samples <- bins.res$count.first.samples
                  for(nn in 1:n.samples){
                     count.samples[indcount,nn] <- count.first.samples[[nn]]
                  }
               }
               colnames(count.samples) <- samples.name
               rownames(count.samples) <- allbins
            }

            ##==========================================
            ## Regularization Path
            ## Use Dichotomie - Use REFIT
            ##==========================================
            loss.weights <- sum(NN)*len/1e9
            param$loss_weights <- as.matrix(loss.weights)

            if(n.samples==1 | mergerefit) { # use the 1 dim fit against the one sample or the pool
               if(!mergerefit & !is.null(samples.name)){
                  count.mono <- count.samples # the count of one given sample
               } else {
                  count.mono <- count # the pool
               }
               respath <- regularization_path(graph, count.mono, param, max_isoforms, delta,
                                              fast_guess=1, iterpoisson=2000, tolpoisson=1e-6, verbosepath, verbose)
               path.set=respath$path.set
            } else { # use group-lasso approach
               # 1) collect paths --------------------------------------------------
               # --> on the sum counts 
               # --> on the individual counts
               n.samples.eff <- n.samples+1
               count.samples.eff <- cbind(count, count.samples) 
               all.collpath <- mclapply(1:n.samples.eff, FUN=function(jj){
                                        collpath <- collect_path_grouplasso(graph, count.samples.eff[,jj,drop=F], param, max_isoforms, delta, 
                                                                            fast_guess=1, iterpoisson=2000, tolpoisson=1e-6, verbosepath, verbose)
                                        pm=matrix(unlist(collpath$path.set), nrow=n.nodes)
                                        bm=matrix(unlist(collpath$beta.avantrefit), nrow=1)
                                        indup <- duplicated(t(pm), fromLast=TRUE)
                                        pm <- pm[,!indup, drop=F]
                                        bm <- bm[,!indup, drop=F]
                                        return(list(pm, bm))
                                       }, mc.cores=mc.cores)
               path.coll <- matrix(unlist(lapply(all.collpath, FUN=function(gg) gg[[1]])), nrow=n.nodes)
               beta.coll <- matrix(unlist(lapply(all.collpath, FUN=function(gg) gg[[2]])), nrow=1)
               indup <- duplicated(t(path.coll), fromLast=TRUE)
               path.coll <- path.coll[,!indup, drop=F]
               beta.coll <- beta.coll[,!indup, drop=F]
               # je peux sans doute faire mieux pour l'initialisation -- 
               # prendre les beta de chaque sample s'il existe le chemin ......
               respath <- regularization_path_grouplasso_dich(path.coll=path.coll, beta.coll=beta.coll,
                                                              count.samples=count.samples, loss.weights=loss.weights,
                                                              delta=delta, iterpoisson=2000, tolpoisson=1e-6, tolpoisson_refit=1e-20,
                                                              REFIT=TRUE, verbosepath=verbosepath, verbose=verbose,
                                                              n.samples=n.samples, NN=NN, len=len, max_isoforms=max_isoforms, mc.cores=mc.cores)
               path.tmp=path.coll
            }
            beta.refit=respath$beta.refit # ok in FPKM with individual normalization
            size.set=respath$size.set
            loss.set=respath$loss.set
            #lambda.set=respath$lambda.set
            beta.avantrefit=respath$beta.avantrefit

            if(verbose==1) print('MODEL SELECTION ...')
            ##================
            ## Model Selection:
            ## BIC Criterion
            ##================			
            select <- modelselection(BICcst=BICcst, loss.set=loss.set,
                                     beta.refit=beta.refit, beta.avantrefit=beta.avantrefit,
                                     size.set=size.set, n.nodes=n.nodes,
                                     n.samples=n.samples, mergerefit=mergerefit)

            # beta in FPKM and path
            beta.select <- beta.refit[[select]] # abundance values in FPKM
            if(n.samples==1 | mergerefit) {
               path.tmp <- path.set[[select]]
               if(mergerefit) {  # AMELIORER LES CONDITIONS ?
                  beta.select <- refitSamples(count.samples, path.tmp, beta.select, NN, len, delta, 
                                              iterpoisson=2000, tolpoisson=1e-20, verbosepath, mc.cores=mc.cores) # multi-samples refit
               } else {
                  beta.select <- t(beta.select)
               }
            }
            if(ncol(path.tmp)==0) path.tmp <- matrix(0, n.nodes, 1) # 2015-03-02
            rownames(beta.select) <- samples.name
            path.select <- apply(path.tmp, 2, FUN=function(pp){
                                 onep <- rep(0, n.exons)
                                 ind.exons <- unique(unlist(strsplit(allbins[which(pp!=0)], split='-')))
                                 ind.exons <- sort(as.numeric(ind.exons))
                                 onep[ind.exons] <- 1
                                 return(onep)})
            if(!is.matrix(path.select)) path.select <- matrix(path.select, nrow=n.nodes)  # (for when only 1 bins) 
            if(ncol(path.select)==0) path.select <- matrix(0, nrow=n.nodes) # 2015-03-02
            # 2015-03-02
            # no need to check for repeating paths for paired anymore
            # cutoff is now in the writing part
            npaths <- ncol(path.select)
            ind.exons.all <- lapply(1:npaths, FUN=function(tt){ which(path.select[,tt]!=0)})
            beta.raw <- sapply(1:npaths, FUN=function(qq){ # beta in RAW COUNTS
                               ind.exons <- ind.exons.all[[qq]]
                               if(eff.paired){
                                  return(pmax(0,(sum(len.exons[ind.exons]) - frag + 1)*beta.select[,qq]*NN/1e9)) # expected raw counts for paired-end
                               }
                               if(!eff.paired){
                                  return(pmax(0,(sum(len.exons[ind.exons]) - readlen + 1)*beta.select[,qq]*NN/1e9)) # expected raw counts for single-end
                                 }})
            if(!is.matrix(beta.raw)) beta.raw <- matrix(beta.raw, ncol=npaths)
            rownames(beta.raw) <- rownames(beta.select)

            if(length(beta.select)>0){ # UTILE ? 
               if(verbose==1) print('WRITE in OUTPUT ...')
               ##================
               ## WRITE OUTPUT
               ##================
               WriteOutput(tophat.exons=tophat.exons, npaths=npaths, ind.exons.all=ind.exons.all,
                           chr=chr, nom=nom, strand=strand,
                           beta.select=beta.select, beta.raw=beta.raw, cutoff=cutoff,
                           n.samples=n.samples, samples.name=samples.name,  
                           output.type=output.type, outf=outf, outf.fpkm=outf.fpkm, outf.count=outf.count)
            }
         }
         timer <- toc()
         beta.fpkm[[mm]] <- beta.select
         beta.expcount[[mm]] <- beta.raw
         timer.all[mm] <- timer
         strand.gr <- strand
         if(strand.gr=='.') strand.gr='*'
         gr <- GRanges(chr,
                       IRanges(tophat.exons[,1], tophat.exons[,2], names=rep(nom, n.exons)),
                       strand=strand.gr,
                       read.count=count.exons,
                       transcript=path.select)
         granges[[mm]] <- gr
      }

      close(inpf)
      ccl <- lapply(outf, close)
      if(!is.null(samples.name) & output.type=='table') {
         close(outf.fpkm)
         close(outf.count)
      }

      return(list(transcripts=granges, abundancesFPKM=beta.fpkm, expected.counts=beta.expcount, timer=timer.all))

}
