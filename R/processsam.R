processsam <- function(sam.file, prefix, annot='', minReadNum=10, minCvgCut=0.25, verbose=0)
{
   .C('ffProcesssam',
      as.character(sam.file),
      as.character(prefix),
      as.character(annot),
      as.character(minReadNum),
      as.character(minCvgCut),
      as.character(verbose),
      PACKAGE='flipflop')
}
