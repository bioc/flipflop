write.transcript <- function(output.file, chr, name, 
                             strand, transcript.start, transcript.end, 
                             RPKM, num){
  cat(
      chr,
      'Flipflop',
      'transcript',
      transcript.start, transcript.end, 
      '.', strand, '.',
      paste('gene_id', paste("\"", name, "\"", ";", sep=""),
	'transcript_id', paste("\"", paste(name, num, sep="."),"\"", ";", sep=""),
	'FPKM', paste("\"", RPKM, "\"", ";", sep="")
	),
      sep='\t', append=TRUE, '\n', file=output.file)
}

write.exons <- function(output.file, chr, name, 
                        strand, pos.start, pos.stop, 
                        num, number, RPKM){
  cat(
      chr,
      'Flipflop',
      'exon',
      pos.start, pos.stop,
      '.', strand, '.',
      paste('gene_id', paste("\"", name, "\"", ";", sep=""),
	'transcript_id', paste("\"", paste(name, num, sep="."), "\"", ";", sep=""),
	'exon_number', paste("\"", number, "\"", ";", sep=""),
	'FPKM', paste("\"", RPKM, "\"", ";", sep="")
	),
      sep='\t', append=TRUE, '\n', file=output.file)
}

