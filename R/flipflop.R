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
##   \item{data.file}{[input] Input alignment file in SAM format. The SAM file must be sorted according to chromosome name and starting position. If you use a multi-sample strategy, please give a single SAM file that results from the "samtools merge" command with the "-r" option (i.e attach RG tag infered from file name).}
##   \item{out.file}{[output] Output gtf file storing the structure of the transcripts which are found to be expressed together with their abundances (in FPKM and expected count).}
##   \item{output.type}{[output] Type of output when using several samples simultaneously. When equal to "gtf" the output corresponds to a gtf file per sample with specific abundances. When equal to "table" the output corresponds to a single gtf file storing the structure of the transcripts, and an associated table with the abundances of all samples (transcripts in rows and samples in columns) Default "gtf".}
##   \item{annot.file}{[input optional] Optional annotation file in BED12 format. If given, exon boundaries will be taken from the annotations. The BED file should be sorted according to chromosome name and starting position of transcripts.}
##   \item{samples.file}{[multi-samples] Optional samples file (one line per sample name). The names should be the one present in the RG tag of the provided SAM file.}
##   \item{mergerefit}{[multi-samples] If TRUE use a simple refit strategy instead of the group-lasso. Default FALSE.}
##   \item{paired}{[paired] Boolean for paired-end reads. If FALSE your reads will be considered as single-end reads. Default FALSE.}
##   \item{frag}{[paired] Mean fragment size. Only used if paired is set to TRUE. Default 400.}
##   \item{std}{[paired] Standard deviation of fragment size. Only used if paired is set to TRUE. Default 20.}
##   \item{OnlyPreprocess}{[pre-processing] Boolean for performing only the pre-processing step. Output is two files: one file '.instance' and one other file '.totalnumread'. Default FALSE.}
##   \item{preprocess.instance}{[pre-processing] Give directly the pre-processed '.instance' input file created when using the OnlyPreprocess option. If non empty, the data.file and annot.file fields are ignored.}
##   \item{minReadNum}{[pre-processing] The minimum number of clustered reads to output. Default 40. If you give an annotation file it will be the minimum number of mapped reads to process a gene.}
##   \item{minFragNum}{[pre-processing] The minimum number of mapped read pairs to process a gene. Only used if paired is TRUE. Default 20.}
##   \item{minCvgCut}{[pre-processing] The fraction for coverage cutoff, should be between 0-1. A higher value will be more sensitive to coverage discrepancies in one gene. Default 0.05.}
##   \item{minJuncCount}{[pre-processing] The minimum number of reads to consider a junction as valid. Default 1.}
##   \item{NN}{[pre-processing] Total number of mapped fragments. Optional. If given the number of mapped fragments is not read from the '.totalnumread' file.}
##   \item{sliceCount}{[parallelization] Number of slices in which the pre-processing '.instance' file is cut. It creates several instance files with the extension '_jj.instance' where jj is the number of the slice. If you set OnlyPreprocess to TRUE, it will create those slices and you can run FlipFlop independently afterwards on each one of the slice. Default 1.}
##   \item{mc.cores}{[parallelization] Number of cores. If you give sliceCount>1 with OnlyPreprocess set to FALSE, it will distribute the slices on several cores. If you give a preprocess.instance file as input (which might be a slice of an original instance file), it will use several cores when you are using a multi-samples strategy. Default 1.}
##   \item{verbose}{[verbosity] Verbosity. Default 0 (little verbosity). Put 1 for more verbosity.}
##   \item{verbosepath}{[verbosity] Verbosity of the optimization part. Default 0 (little verbosity). Put 1 for more verbosity. Mainly used for de-bugging.}
##   \item{cutoff}{[parameter] Do not report isoforms whose expression level is less than cutoff percent of the most expressed transcripts. Default 1.}
##   \item{BICcst}{[parameter] Constant used for model selection with the BIC criterion. Default 50.}
##   \item{delta}{[parameter] Loss parameter, Poisson offset. Default 1e-7.}
##   \item{use_TSSPAS}{Do we restrict the candidate TSS and PAS sites. 1 is yes and 0 is no. Default 0 i.e each exon can possibly starts or ends an isoform.}
##   \item{max_isoforms}{Maximum number of isoforms given during regularization path. Default 10.}
## }
##
## \value{
##  A @list with the following elements:
##  \item{transcripts}{A list storing the structure of the expressed isoforms. The list is a GenomicRangesList object from the GenomicRanges package. Rows correspond to exons. On the left hand side each exon is described by the gene name, the chromosome, its genomic position on the chromosome and the strand. Transcripts are described on the right hand side. Every transcript is a binary vector where an exon is labelled by 1 if it is included in the transcript.}
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
                     output.type='gtf',
                     annot.file='',
                     samples.file='',
                     mergerefit=FALSE,
                     paired=FALSE,
                     frag=400,
                     std=20,
                     OnlyPreprocess=FALSE,
                     preprocess.instance='',
                     minReadNum=40,
                     minFragNum=20,
                     minCvgCut=0.05, # 2015-01-22
                     minJuncCount=1,
                     sliceCount=1,
                     mc.cores=1, # 2015-03-16
                     NN='',
                     verbose=0,
                     verbosepath=0,
                     cutoff=1, # 2015-01-22
                     BICcst=50,
                     delta=1e-7, # 2015-01-22
                     use_TSSPAS=0,
                     max_isoforms=10
                     ){

   options(scipen=20) 
   
   ##================================================================================
   ## Process alignment data using functions of IsoLasso (Li et al. 2011)
   ##================================================================================
   ## Pre-Process SAM input file: 
   if(preprocess.instance==''){
      print('PRE-PROCESSING sam file ....')
      data.file <- path.expand(path=data.file) # In case there is a '~' in the input path
      if(data.file==''){ print('Did you forget to give a SAM file?') ; return(NULL) } 
      annot.file <- path.expand(path=annot.file)
      samples.file <- path.expand(path=samples.file)
      bn <- basename(data.file)
      prefix <- sub('[.][^.]*$', '', bn)
      processsam(data.file, prefix, annot=annot.file, samples=samples.file,
                 paired=paired, minReadNum=minReadNum, minCvgCut=minCvgCut, minJuncCount=minJuncCount,
                 verbose=verbose, sliceCount=sliceCount) # Pre-Processing step (mainly from IsoLasso software)
      print('DONE !')
   }

   ## Continue the job:
   if(!OnlyPreprocess){
      
      # Total number of mapped reads, if not provided
      if(NN==''){ # if NN (total number of mapped fragments) is not given, read it in the file created previously
         if(preprocess.instance==''){
            numrfn <- paste(prefix, '.totalnumread', sep="")
            numrf <- file(numrfn,'r') # Read the file where the total number of reads is stored
         }
         if(preprocess.instance!=''){
            go1 <- dirname(preprocess.instance)
            go2 <- sub('[.][^.]*$', '',basename(preprocess.instance))
            goname <- paste(go2,'.totalnumread',sep="")
            if( length(list.files(path=go1, pattern=goname))>0 ){ # you find the associated .totalnumread file
               numrfn <- paste(go1,goname,sep="/")
            } else { # you look if the reason why you did not find it is because you are using a splitted instance file 
               go3 <- sub('(_)[0-9]+(.instance)', '', basename(preprocess.instance), perl=TRUE)
               gonamebis <- paste(go3,'.totalnumread',sep="")
               numrfn <- paste(go1,gonamebis,sep="/")
            }
            numrf <- file(numrfn,'r')
         }
         nn <- scan(numrf, what=character(0), nlines=1, quiet=TRUE, skip=1) # Total number of mapped reads
         if(paired==FALSE){
            NN <- as.integer(nn[1])
         }
         if(paired==TRUE){
            npair <- scan(numrf, what=character(0), nlines=1, quiet=TRUE, skip=1)
            Npair <- as.integer(npair[1])
            NN <- Npair # Total number of mapped fragments
         }
         close(numrf)
      }

      # Samples file
      samples.name <- NULL
      n.samples <- 1
      if(samples.file!=''){
         samples <- read.table(samples.file)
         samples.name <- sort(unique(samples[,1]))
         n.samples <- length(samples.name)
         print(paste('The sample names are:', paste(samples.name, collapse=" ")))
      }

      # Output file & Update NN
      outf.fpkm <- NULL
      outf.count <- NULL
      ind.samples <- NULL
      if(is.null(samples.name)){
         outf <- list(file(out.file, 'w'))
      } else {
         go1 <- dirname(out.file)
         go2 <- sub('[.][^.]*$', '', basename(out.file))
         go <- paste(go1,go2,sep='/')
         tmp.names.all <- sapply(strsplit(nn[-1], split='='), FUN=function(k) k[1]) # all sample names in the instance file
         ind.samples <- sapply(samples.name, FUN=function(ii) which(tmp.names.all==ii)) # keep the ones given in the sample file
         if(output.type=='table') {
            outf <- list(file(out.file, 'w'))
            outf.fpkm <- file(paste(go, 'fpkm', sep='_'), 'w')
            outf.count <- file(paste(go, 'count', sep='_'), 'w')
            cat("Name", as.vector(samples.name), sep="\t", '\n', file=outf.fpkm)
            cat("Name", as.vector(samples.name), sep="\t", '\n', file=outf.count)
         } else {
            outf <- lapply( samples.name, FUN=function(ss) file(paste(go,'_',ss,'.gtf',sep=""),'w') )
         }
         if(!paired) {
            nn.parse <- strsplit(nn[-1], split='=')
            NN <- sapply(nn.parse[ind.samples], FUN=function(k) as.integer(k[2]))
         } else {
            npair.parse <- strsplit(npair[-1], split='=')
            NN <- sapply(npair.parse[ind.samples], FUN=function(k) as.integer(k[2]))
         }
      }

      # Use given preprocess.instance input file:
     doparal <- FALSE 
      if(preprocess.instance!='') {
         inpf <- file(preprocess.instance, 'r')
      } else {
         infn <- paste(prefix, '.instance', sep='')
         # Or read the preprocess.instance file just created:
         if(sliceCount<=1) { # you read the global instance file
            inpf <- file(infn, 'r')
         } else { # you will read the slices of the instance file
            slicefn <- paste(prefix,'.slicecount',sep="")
            slicef <- file(slicefn,'r')
            numsl <- scan(slicef, what=character(0), nlines=1, quiet=TRUE, skip=1) # number of effective slices
            sliceCount.eff <- as.integer(numsl[1])
            close(slicef)
            doparal <- TRUE
         }
      }

      if(!doparal) {
         res.flipflop <- WrapperFlipFlop(inpf=inpf, # the input open file
                                         outf=outf, # the output open file
                                         outf.fpkm=outf.fpkm,
                                         outf.count=outf.count,
                                         output.type=output.type,
                                         NN=NN,
                                         samples.name=samples.name,
                                         n.samples=n.samples,
                                         ind.samples=ind.samples,
                                         mergerefit=mergerefit,
                                         paired=paired,
                                         frag=frag,
                                         std=std,
                                         minReadNum=minReadNum,
                                         minFragNum=minFragNum,
                                         mc.cores=mc.cores,
                                         verbose=verbose,
                                         verbosepath=verbosepath,
                                         cutoff=cutoff,
                                         BICcst=BICcst,
                                         delta=delta,
                                         use_TSSPAS=use_TSSPAS,
                                         max_isoforms=max_isoforms
                                         )
      } else {
         mc.cores.inside <- 1
         if( n.samples > 1 & (mc.cores - sliceCount.eff >= 2*n.samples) ) { mc.cores.inside <- mc.cores - sliceCount.eff} # you have at least 2 cores per slice available
         res <- mclapply(1:sliceCount.eff, FUN=function(ss){
                         infn.slice <- paste(prefix,'_',ss,'.instance', sep='')
                         inpf.slice <- file(infn.slice, 'r')
                         outf.fpkm.slice <- NULL
                         outf.count.slice <- NULL
                         if(!is.null(outf.fpkm)){
                            go1 <- dirname(out.file)
                            go2 <- sub('[.][^.]*$', '', basename(out.file))
                            go <- paste(go1,go2,sep='/')
                            outf.fpkm.slice <- file(paste(go, 'fpkm', ss, sep='_'), 'w')
                            outf.count.slice <- file(paste(go, 'count', ss, sep='_'), 'w')
                         }
                         outf.slice <- lapply(outf, FUN=function(oo) {
                                              go <- sub('[.][^.]*$','',summary(oo)$description)
                                              return(file(paste(go, '_', ss,'.gtf', sep=''), 'w')) })
                         res.flipflop <- WrapperFlipFlop(inpf=inpf.slice, # the input open file
                                                         outf=outf.slice, # the output open file
                                                         outf.fpkm=outf.fpkm.slice,
                                                         outf.count=outf.count.slice,
                                                         output.type=output.type,
                                                         NN=NN,
                                                         samples.name=samples.name,
                                                         n.samples=n.samples,
                                                         ind.samples=ind.samples,
                                                         mergerefit=mergerefit,
                                                         paired=paired,
                                                         frag=frag,
                                                         std=std,
                                                         minReadNum=minReadNum,
                                                         minFragNum=minFragNum,
                                                         mc.cores=mc.cores.inside,
                                                         verbose=verbose,
                                                         verbosepath=verbosepath,
                                                         cutoff=cutoff,
                                                         BICcst=BICcst,
                                                         delta=delta,
                                                         use_TSSPAS=use_TSSPAS,
                                                         max_isoforms=max_isoforms
                                                         )
                         #unlink(infn.slice)  # provoque un bug je ne sais pas pourquoi... !?
                         }, mc.cores=mc.cores)
         # concatenate output in a single one and erase the slices
         if(!is.null(outf.fpkm)){
            go1 <- dirname(out.file)
            go2 <- sub('[.][^.]*$', '', basename(out.file))
            go <- paste(go1,go2,sep='/')
            system( paste('cat', paste(go, 'fpkm', 1:sliceCount.eff, sep="_", collapse=" "), '>', paste(go,'fpkm',sep="_"), sep=" ") )
            system( paste('cat', paste(go, 'count', 1:sliceCount.eff, sep="_", collapse=" "), '>', paste(go,'count',sep="_"), sep=" ") )
            lapply(1:sliceCount.eff, FUN=function(ss){ unlink(paste(go, 'fpkm', ss, sep="_")) ; unlink(paste(go, 'count', ss, sep="_"))})
         }
         aa <- lapply(outf, FUN=function(oo) {
                      go <- sub('[.][^.]*$','',summary(oo)$description)
                      system( paste('cat', paste(go, '_', 1:sliceCount.eff,'.gtf', sep="", collapse=" "), '>', paste(go,'.gtf',sep=""), sep=" ") ) 
                      lapply(1:sliceCount.eff, FUN=function(ss){ unlink(paste(go, '_', ss, '.gtf', sep=""))})
                                         })
         res.flipflop <- unlist(res, recursive=FALSE)
      }
      granges <- res.flipflop$transcripts
      beta.fpkm <- res.flipflop$abundancesFPKM
      beta.expcount <- res.flipflop$expected.countd
      timer.all <- res.flipflop$timer
   }
   ## Delete processed sam file
   if(preprocess.instance=='' & !OnlyPreprocess){
      unlink(numrfn)
      unlink(infn)
      if(doparal){ 
         unlink(slicefn)
         lapply(1:sliceCount.eff, FUN=function(ss) unlink(paste(prefix,'_',ss,'.instance', sep="")))
      }
   }
   if(OnlyPreprocess==FALSE){
      return(list(transcripts=granges, abundancesFPKM=beta.fpkm, expected.counts=beta.expcount, timer=timer.all))
   }
}
