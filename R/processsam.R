processsam <- function(sam.file, prefix, annot='', samples='', paired=FALSE, minReadNum=10, minCvgCut=0.25, minJuncCount=1, verbose=0)
{
   .C('ffProcesssam',
      as.character(sam.file),
      as.character(prefix),
      as.character(annot),
      as.character(samples), # 2015-01-15
      as.character(paired),
      as.character(minReadNum),
      as.character(minCvgCut),
      as.character(minJuncCount),
      as.character(verbose),
      PACKAGE='flipflop')
}
