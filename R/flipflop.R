# Copyright 2013 Elsa Bernard, Laurent Jacob, Julien Mairal and Jean-Philippe Vert

## This file is part of flipflop.

## flipflop is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## flipflop is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with flipflop.  If not, see <http://www.gnu.org/licenses/>.

#########################################################################/**
## @RdocFunction flipflop
##
## @title "Estimate isoform compositions and abundances"
##
## \description{ This function takes count data (RNA-seq alignment in SAM format) for a given gene as
##  input and estimates which isoforms of the gene are most likely to
##  have generated this set of counts. It is based on a Poisson
##  likelihood penalized by an l1 norm as explained in Bernard et al.,
##  2013.}
##
## @synopsis
##
## \arguments{
##   \item{data.file}{Input alignment file in SAM format. The SAM file must be sorted according to chromosome name and starting position.}
##   \item{out.file}{Output gtf file storing the structure of the transcripts which are found to be expressed together with their abundances (in FPKM and expected count).}
##   \item{annot.file}{Optional annotation file in BED12 format. If given, exon boundaries will be taken from the annotations. The BED file should be sorted according to chromosome name and starting position of transcripts.}
##   \item{paired}{Boolean for paired-end reads. If FALSE your reads will be considered as single-end reads. Default FALSE.}
##   \item{frag}{Mean fragment size. Only used if paired is set to TRUE. Default 400.}
##   \item{std}{Standard deviation of fragment size. Only used if paired is set to TRUE. Default 20.}
##   \item{minReadNum}{[Pre-processing] The minimum number of clustered reads to output. Default 40. If you give an annotation file it will be the minimum number of mapped reads to process a gene.}
##   \item{minFragNum}{[Pre-processing] The minimum number of mapped read pairs to process a gene. Only used if paired is TRUE. Default 20.}
##   \item{minCvgCut}{[Pre-processing] The fraction for coverage cutoff, should be between 0-1. A higher value will be more sensitive to coverage discrepancies in one gene. Default 0.25.}
##   \item{minJuncCount}{[Pre-processing] The minimum number of reads to consider a junction as valid. Default 1.}
##   \item{verbose}{Verbosity. Default 0 (little verbosity). Put 1 for more verbosity.}
##   \item{verbosepath}{Verbosity of the optimization part. Default 0 (little verbosity). Put 1 for more verbosity.}
##   \item{max_isoforms}{Maximum number of isoforms given during regularization path. Default 10.}
##   \item{use_TSSPAS}{Do we restrict the candidate TSS and PAS sites. 1 is yes and 0 is no. Default 0 i.e each exon can possibly starts or ends an isoform.}
##   \item{cutoff}{For paired-end reads do not report isoforms whose expression level is less than cutoff percent of the most expressed transcripts. Not active if paired is FALSE. Default 5.}
##   \item{BICcst}{Constant used for model selection with the BIC criterion. Default 50.}
##   \item{OnlyPreprocess}{Boolean for performing only the pre-processing step. Output is two files: one file '.instance' and one other file '.totalnumread'. Default FALSE.}
##   \item{preprocess.instance}{Give directly the pre-processed '.instance' input file created when using the OnlyPreprocess option. If non empty, the data.file and annot.file fields are ignored.}
##   \item{NN}{Total number of mapped fragments. Optional. If given the number of mapped fragments is not read from the '.totalnumread' file.}
## }
##
## \value{
##  A @list with the following elements:
##  \item{transcripts}{A list storing the structure of the expressed isoforms. The list is a GRangesList object from the GenomicRanges package. Rows correspond to exons. On the left hand side each exon is described by the gene name, the chromosome, its genomic position on the chromosome and the strand. Transcripts are described on the right hand side. Every transcript is a binary vector where an exon is labelled by 1 if it is included in the transcript.}
##  \item{abundancesFPKM}{A list storing the abundances of the expressed isoforms in FPKM unit. Each element of the list is a vector whose length is the number of expressed transcripts listed in the above 'transcripts' object.}
##  \item{expected.counts}{A similar list as 'abundancesFPKM' but storing the expected fragment counts for each expressed isoforms.}
##  \item{timer}{A vector with the computation time in seconds for each gene.}
## }
##
## @author
##
## @examples "../inst/extdata/flipflop.Rex"
##
##*/########################################################################

flipflop <- function(data.file, 
                     out.file='FlipFlop_output.gtf',
                     annot.file='',
                     paired=FALSE,
                     frag=400,
                     std=20,
                     minReadNum=40,
                     minFragNum=20,
                     minCvgCut=0.25,
                     minJuncCount=1,
                     verbose=0,
                     verbosepath=0,
                     max_isoforms=10,
                     use_TSSPAS=0,
                     cutoff=5,
                     BICcst=50,
                     OnlyPreprocess=FALSE,
                     preprocess.instance='',
                     NN=''
                     ){

   options(scipen=20) 
   ##================================================================================
   ## Process alignment data using functions of IsoLasso (Li et al. 2011)
   ##================================================================================
   # Pre-Process sam input file: 
   if(preprocess.instance=='' | OnlyPreprocess==TRUE){
      print('PRE-PROCESSING sam file ....')
      data.file <- path.expand(path=data.file) # In case there is a '~' in the input path 
      annot.file <- path.expand(path=annot.file)
      bn <- basename(data.file)
      prefix <- sub('[.][^.]*$', '', bn)
      processsam(data.file, prefix, annot=annot.file, paired=paired, minReadNum=minReadNum, minCvgCut=minCvgCut, minJuncCount=minJuncCount, verbose=verbose) # Pre-Processing step (mainly from IsoLasso software)
      print('DONE !')
   }

   # Continue the job:
   if(OnlyPreprocess==FALSE){
      # Read the preprocess.instance file just created:
      if(preprocess.instance==''){
         infn <- paste(prefix, '.instance', sep='')
         inpf <- file(infn, 'r')
      }
      # Or use given preprocess.instance input file. 
      if(preprocess.instance!=''){
         inpf <- file(preprocess.instance, 'r')
      }
      # Total number of mapped reads, if not provided
      if(NN==''){ # if NN (total number of mapped fragments) is not given, read it in the file created previously
         if(preprocess.instance==''){
            numrfn <- paste(prefix, '.totalnumread', sep='')
            numrf <- file(numrfn,'r') # Read the file where the total number of reads is stored
         }
         if(preprocess.instance!=''){
            go1 <- dirname(preprocess.instance)
            go2 <- sub('[.][^.]*$', '',basename(preprocess.instance))
            numrfn <- paste(go1,'/',go2,'.totalnumread',sep="")
            numrf <- file(numrfn,'r')
         }
         if(paired==FALSE){
            NN <- scan(numrf, what=integer(0), nlines=1, quiet=TRUE, skip=1) # Total number of mapped reads
         }
         if(paired==TRUE){
            Nall <- scan(numrf, what=integer(0), nlines=1, quiet=TRUE, skip=1)
            Npair <- scan(numrf, what=integer(0), nlines=1, quiet=TRUE, skip=1)
            NN <- Npair # Total number of mapped fragments
         }
         close(numrf)
      }
      outf <- file(out.file, 'w')

      ##=======================================
      ## Set up parameters
      ##=======================================
      delta <- 1e-5 # Smoothing parameter for the Poisson Model
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
      #granges <- GRangesList() 
      granges <- list() # GRangesList does not work because the number of metadata columns varies
      mm <- 0
      scan(inpf, what=character(0), nlines=3, quiet=TRUE)

      while(1){
         mm <- mm+1

         ## Name-Instance
         name <- scan(inpf, what=character(0), nlines=1, quiet=TRUE)
         print(name)
         nom <- paste('Inst', name[2], sep='')

         if(length(name)==0){
            break
         }

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
            listref <- scan(inpf, what=character(0), nlines=nref, quiet=TRUE) # so far we don't use reference information for more than exon boundaries
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
            readtype <- scan(inpf, what=numeric(0), nlines=ntype, quiet=TRUE)
            dim(readtype) <- c((n.exons+2), ntype)
            readtype <- t(readtype)
         }
         if(ntype==0){
            readtype <- matrix(0,1,(n.exons+1))
         }

         ## Paired-End type
         petype <- scan(inpf, what=character(0), nlines=1, quiet=TRUE)
         npereads <- as.integer(petype[2])
         npetype <- as.integer(petype[3])
         # Skip instance when not enough mapped reads or mapped fragments (but don't skip if your reads are supposed to be paired but there is more than 1000 single-end reads):
         if(nreads<minReadNum | (paired==TRUE & npereads<minFragNum & nreads<1000)){
            if(npetype>0){
               scan(inpf, what=character(0), nlines=(2*npetype), quiet=TRUE) # Skip paired-end type information
            }
            next
         }
         # Consider reads as single-end (npetype==0 means there is no pair and for npetype>300 the paired graph creation is too long):
         eff.paired <- FALSE
         if(npetype>0){
            if(paired==FALSE | npetype>=300) {
               scan(inpf, what=character(0), nlines=(2*npetype), quiet=TRUE) # Skip paired-end type information
            } else {
            # Paired-end reads:
               eff.paired <- TRUE
               # Construct the bins (consider pairs as 'long single' with some heuristics if necessary)
               bins.res <- bins_paired(inpf=inpf, npetype=npetype, tophat.exons=tophat.exons, count.exons=count.exons, 
                                       len.exons=len.exons, n.exons=n.exons, readtype=readtype, gap=gap, std=std) 
               keeplolo <- bins.res$longtype
               count.first <- bins.res$count.first
               VALID <- (length(keeplolo)>0)
            }
         }

         tic()
         ## MULTI-EXONS
         if(n.exons>1){
            ##==========================================
            ## Create Splicing Graph
            ## NODE:=READTYPE:=BIN
            ##==========================================
            if(!eff.paired){
               #lolo <- lapply(1:nrow(readtype), FUN=function(v) which(readtype[v,1:n.exons]==1))
               ### sanity check 1 --- duplicated lolo
               ## -------------------------------------
               #list.dup <- lolo[duplicated(lolo)]
               #to.remove <- c()
               #for(dd in unique(list.dup)){
               #   # indice des dupliques
               #   ind.dup <- which(sapply(lolo, FUN=function(k) identical(k,dd)))
               #   readtype[ind.dup[1],n.exons+1] <- sum(readtype[ind.dup,n.exons+1])
               #   to.remove <- c(to.remove, ind.dup[2:length(ind.dup)])
               #}
               #if(!is.null(to.remove)){
               #   readtype <- readtype[-to.remove,,drop=F]
               #}
               lolo <- lapply(1:nrow(readtype), FUN=function(v) which(readtype[v,1:n.exons]==1))
               ## sanity check 2 --- valid bins
               # -----------------------------------
               validbins <- sapply(lolo, FUN=function(v){
                                   if(length(v)<=2){
                                      return((sum(len.exons[v])>=readlen))
                                   }
                                   if(length(v)>=3){
                                      vv <- v[2:(length(v)-1)]
                                      return(((sum(len.exons[v])>=readlen) & (sum(len.exons[vv])<(readlen-1)))) 
                                   }}) 
               VALID <- (sum(validbins)>0)
            }

            if(nreads == 0 | !VALID) { # don't go if you do not have any bins
               path.select <- matrix(0, nrow=n.exons, ncol=1)
               beta.select <- 0
            } else {
               # Create bins and junctions
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

               ##==========================================
               ## Regularization Path
               ## Use Dichotomie - Use REFIT
               ##==========================================
               loss.weights <- NN*len/1e9
               param$loss_weights <- as.matrix(loss.weights)

               respath <- regularization_path(graph, count, param, max_isoforms, delta, fast_guess=1, iterpoisson=2000, tolpoisson=1e-9, verbosepath, verbose)
               beta.refit=respath$beta.refit
               path.set=respath$path.set
               size.set=respath$size.set
               loss.set=respath$loss.set

               if(verbose==1){
                  print('MODEL SELECTION ...')
               }
               ##================
               ## Model Selection:
               ## BIC Criterion
               ##================			
               beta.select <- 0 # abundance values in FPKM
               beta.raw <- 0 # expected fragment 'raw' counts 
               path.select <- matrix(0, n.exons)
               select <- which.min(loss.set + BICcst*size.set*log(n.nodes))
               trueind <- which(beta.refit[[select]] > 0)
               if(length(trueind)>0){
                  beta.select <- beta.refit[[select]][trueind]
                  path.tmp <- matrix(path.set[[select]][,trueind],nrow=n.nodes)
                  rownames(path.tmp) <- allbins
                  if(eff.paired){ # For paired-end check that there is no repeating paths:
                     if(length(beta.select)==1){
                        ind.exons <- unique(unlist(strsplit(allbins[which(path.tmp[,1]!=0)], split='-')))
                        ind.exons <- sort(as.numeric(ind.exons))
                        path.select <- matrix(0, nrow=n.exons, ncol=1)
                        path.select[ind.exons,1] <- 1
                     }
                     if(length(beta.select)>1){
                        listpath <- lapply(1:ncol(path.tmp), FUN=function(v) sort(as.numeric(unique(unlist(strsplit(allbins[which(path.tmp[,v]!=0)],split='-'))))))
                        tutu <- duplicated(listpath)
                        if(sum(tutu)>0){
                           indrep <- which(tutu==TRUE)
                           for(gg in indrep){
                              koi <- which(sapply(listpath,FUN=function(v) identical(v, listpath[[gg]])))
                              beta.select[koi[koi!=gg]] <- sum(beta.select[koi])
                           }
                           listpath <- listpath[tutu==FALSE]
                           beta.select <- beta.select[tutu==FALSE]
                        }
                        path.select <- sapply(listpath, FUN=function(ll){
                                              onep <- rep(0, n.exons)
                                              onep[ll] <- 1
                                              return(onep)})
                        #-----------------------------#
                        ### CUT-OFF for low value (5% default) 
                        keepend <- which(beta.select > max(beta.select)*cutoff/100)
                        beta.select <- beta.select[keepend]
                        path.select <- path.select[,keepend, drop=F]
                        #-----------------------------#
                     }
                  }
                  if(!eff.paired){
                     path.select <- apply(path.tmp, 2, FUN=function(pp){
                                          onep <- rep(0, n.exons)
                                          ind.exons <- unique(unlist(strsplit(allbins[which(pp!=0)], split='-')))
                                          ind.exons <- sort(as.numeric(ind.exons))
                                          onep[ind.exons] <- 1
                                          return(onep)})
                  }
                  npaths <- length(beta.select)
                  ind.exons.all <- lapply(1:npaths, FUN=function(tt){ which(path.select[,tt]!=0)})
                  beta.raw <- sapply(1:npaths, FUN=function(qq){ 
                                     ind.exons <- ind.exons.all[[qq]]
                                     if(eff.paired){
                                       return(max(0,(sum(len.exons[ind.exons]) - frag + 1)*beta.select[qq]*NN/1e9)) # expected raw counts for paired-end
                                     }
                                     if(!eff.paired){
                                        return(max(0,(sum(len.exons[ind.exons]) - readlen + 1)*beta.select[qq]*NN/1e9)) # expected raw counts for single-end
                                     }})
                  if(verbose==1){
                     print('WRITE in OUTPUT ...')
                  }
                  ##================
                  ## WRITE OUTPUT
                  ##================
                  manage.exons <- tophat.exons
                  manage.exons[length(manage.exons)] <- manage.exons[length(manage.exons)]+1
                  for(ll in 1:npaths){
                     ind.exons <- ind.exons.all[[ll]]
                     transcript.start <- tophat.exons[(ind.exons[1]),1]
                     transcript.end <- tophat.exons[tail(ind.exons,1),2]
                     write.transcript(outf, chr, nom, strand, transcript.start, transcript.end, beta.select[ll], beta.raw[ll], ll) 
                     numex <- 0
                     subexons <- as.vector(manage.exons[ind.exons,])
                     ii1 <- duplicated(subexons)
                     ii2 <- duplicated(subexons, fromLast=TRUE)
                     subexons <- subexons[!(ii1|ii2)]
                     dim(subexons) = c(length(subexons)/2,2)
                     for(kk in 1:nrow(subexons)){
                        numex <- numex+1
                        write.exons(outf, chr, nom, strand, (subexons[kk,1]), (subexons[kk,2]-1), ll, numex, 
                                    beta.select[ll], beta.raw[ll])
                     }
                  }
               }
            }
         }
         ## MONO-EXON
         if(n.exons==1){
            if(paired==TRUE & npetype>0){
               count.exons <- npereads
            }
            select <- which.min(c(loss.ll(count.exons,0,delta), loss.ll(count.exons, count.exons, delta) + BICcst))
            if(select==2){
               beta.select <- count.exons*1e9/(NN*len.exons)
               path.select <- matrix(1)
               beta.raw <- count.exons
            }
            if(select==1){
               beta.select <- 0
               path.select <- matrix(0)
               beta.raw <- 0
            }
            if(verbose==1){
               print('WRITE in OUTPUT ...')
            }
            ##================
            ## WRITE OUTPUT
            ##================  
            if(beta.select > 0){
               transcript.start <- tophat.exons[1,1]
               transcript.end <- tophat.exons[1,2]
               write.transcript(outf, chr, nom, strand, transcript.start, transcript.end, beta.select, beta.raw, 1)
               write.exons(outf, chr, nom, strand, transcript.start, transcript.end, 1, 1, beta.select, beta.raw)
            }
         }

         timer <- toc()
         beta.fpkm[[mm]] <- beta.select
         beta.expcount[[mm]] <- beta.raw
         timer.all[mm] <- timer
         strand.gr <- strand
         if(strand.gr=='.'){
            strand.gr='*'
         }
         gr <- GRanges(chr, 
                       IRanges(tophat.exons[,1], tophat.exons[,2], names=rep(nom, n.exons)),
                       strand=strand.gr,
                       read.count=count.exons,
                       transcript=path.select)
         granges[[mm]] <- gr
      }
      close(inpf)
      close(outf)
   }
   ## Delete processed sam file
   if(preprocess.instance=='' & OnlyPreprocess==FALSE){
      unlink(numrfn)
      unlink(infn)
   }
   if(OnlyPreprocess==FALSE){
      return(list(transcripts=granges, abundancesFPKM=beta.fpkm, expected.counts=beta.expcount, timer=timer.all))
   }
}

