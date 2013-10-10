## Load the library
library(flipflop)

## Alignment data file in SAM format 
data.file <- system.file('extdata/vignette-sam.txt', package='flipflop')

## Run flipflop
ff.res	<- flipflop(data.file=data.file, 
                    out.file='FlipFlop_output_example.gtf')

## Names of the result list returned by flipflop
names(ff.res)

## Structure of the expressed isoforms for the first gene
## Rows correspond to exons, with chromosome, genomic position and strand information for each exon
## The metadata columns correspond to the expressed transcripts 
transcripts <- ff.res$transcripts[[1]]
print(transcripts)

## Abundances in FPKM of the expressed isoforms for the first gene
## The length of the vector corresponds to the number of transcripts listed in the 'transcripts' object.
## Each element of the vector is the estimated abundance of the corresponding transcript.
abundances <- ff.res$abundances[[1]]
print(abundances)
