WriteOutput <- function(tophat.exons, npaths, ind.exons.all,
                        chr, nom, strand, 
                        beta.select, beta.raw, cutoff,
                        n.samples, samples.name,
                        output.type, outf, outf.fpkm, outf.count) {


   manage.exons <- tophat.exons
   manage.exons[length(manage.exons)] <- manage.exons[length(manage.exons)]+1

   for(ll in 1:npaths) {

      ind.exons <- ind.exons.all[[ll]]
      transcript.start <- tophat.exons[(ind.exons[1]),1]
      transcript.end <- tophat.exons[tail(ind.exons,1),2]

      if(!is.null(samples.name) & output.type=='table') { # one unique GTF and TABLES
         if( sum(beta.select[,ll])>0 ){ # expressed in at least one sample
            # write the single GTF
            write.transcript.short(outf[[1]], chr, nom, strand, transcript.start, transcript.end, ll)
            numex <- 0
            subexons <- as.vector(manage.exons[ind.exons,])
            ii1 <- duplicated(subexons)
            ii2 <- duplicated(subexons, fromLast=TRUE)
            subexons <- subexons[!(ii1|ii2)]
            dim(subexons) = c(length(subexons)/2,2)
            for(kk in 1:nrow(subexons)){
               numex <- numex+1
               write.exons.short(outf[[1]], chr, nom, strand, (subexons[kk,1]), (subexons[kk,2]-1), ll, numex)
            }
            # write tables of FPKM and COUNT
            cat(paste(nom,ll,sep='.'), beta.select[,ll], sep='\t', '\n', append=T, file=outf.fpkm)
            cat(paste(nom,ll,sep='.'), beta.raw[,ll], sep='\t', '\n', append=T, file=outf.count)
         }
      } else { # one GTF per sample
         for( ss in 1:n.samples ){
            # 2015-01-19 cutoff for both single and paired
            # 2015-03-11 add fraction in output GTF
            indsup <- which( beta.select[ss,] > max(beta.select[ss,])*cutoff/100)
            total <- sum(beta.select[ss,indsup])
            #if( beta.select[ss,ll] > max(beta.select[ss,])*cutoff/100 ){
            if( ll %in% indsup){
               frac <- beta.select[ss,ll]/total
               write.transcript(outf[[ss]], chr, nom, strand, transcript.start, transcript.end, beta.select[ss,ll], beta.raw[ss,ll], frac, ll) 
               numex <- 0
               subexons <- as.vector(manage.exons[ind.exons,])
               ii1 <- duplicated(subexons)
               ii2 <- duplicated(subexons, fromLast=TRUE)
               subexons <- subexons[!(ii1|ii2)]
               dim(subexons) = c(length(subexons)/2,2)
               for(kk in 1:nrow(subexons)){
                  numex <- numex+1
                  write.exons(outf[[ss]], chr, nom, strand, (subexons[kk,1]), (subexons[kk,2]-1), ll, numex, 
                              beta.select[ss,ll], beta.raw[ss,ll], frac)
               }
            }
         }
      }
   }
}


write.transcript <- function(output.file, chr, name, 
                             strand, transcript.start, transcript.end, 
                             FPKM, EXP_COUNT, frac, num){
   cat(
       chr,
       'Flipflop',
       'transcript',
       transcript.start, transcript.end, 
       '.', strand, '.',
       paste('gene_id', paste("\"", name, "\"", ";", sep=""),
             'transcript_id', paste("\"", paste(name, num, sep="."),"\"", ";", sep=""),
             'FPKM', paste("\"", FPKM, "\"", ";", sep=""),
             'EXP_COUNT', paste("\"", EXP_COUNT, "\"", ";", sep=""),
             'frac', paste("\"", frac, "\"", ";", sep="")
             ),
       sep='\t', append=TRUE, '\n', file=output.file)
}

write.exons <- function(output.file, chr, name, 
                        strand, pos.start, pos.stop, 
                        num, number, FPKM, EXP_COUNT, frac){
   cat(
       chr,
       'Flipflop',
       'exon',
       pos.start, pos.stop,
       '.', strand, '.',
       paste('gene_id', paste("\"", name, "\"", ";", sep=""),
             'transcript_id', paste("\"", paste(name, num, sep="."), "\"", ";", sep=""),
             'exon_number', paste("\"", number, "\"", ";", sep=""),
             'FPKM', paste("\"", FPKM, "\"", ";", sep=""),
             'EXP_COUNT', paste("\"", EXP_COUNT, "\"", ";", sep=""),
             'frac', paste("\"", frac, "\"", ";", sep="")
             ),
       sep='\t', append=TRUE, '\n', file=output.file)
}

write.transcript.short <- function(output.file, chr, name, 
                                   strand, transcript.start, transcript.end, 
                                   num){
   cat(
       chr,
       'Flipflop',
       'transcript',
       transcript.start, transcript.end, 
       '.', strand, '.',
       paste('gene_id', paste("\"", name, "\"", ";", sep=""),
             'transcript_id', paste("\"", paste(name, num, sep="."),"\"", ";", sep="")
             ),
       sep='\t', append=TRUE, '\n', file=output.file)
}

write.exons.short <- function(output.file, chr, name, 
                              strand, pos.start, pos.stop, 
                              num, number){
   cat(
       chr,
       'Flipflop',
       'exon',
       pos.start, pos.stop,
       '.', strand, '.',
       paste('gene_id', paste("\"", name, "\"", ";", sep=""),
             'transcript_id', paste("\"", paste(name, num, sep="."), "\"", ";", sep=""),
             'exon_number', paste("\"", number, "\"", ";", sep="")
             ),
       sep='\t', append=TRUE, '\n', file=output.file)
}
