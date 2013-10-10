library(R.utils)
library(flipflop)

## Set default author
author <- "Elsa Bernard, Laurent Jacob, Julien Mairal and Jean-Philippe Vert"
     
setwd("R")
sourceDirectory(".")
rdocFiles <- list.files(".", full.names=TRUE)

## Compile the Rdoc:s into Rd files (saved in the destPath directory)
destPath <- "../man"
destPath <- Arguments$getWritablePath(destPath)
dummy <- lapply(rdocFiles, Rdoc$compile, destPath=destPath)

## List the generated Rd files
rdFiles <- list.files(destPath, full.names=TRUE)
print(rdFiles)
     
## Show one of the files
## file.show(rdFiles[1])
     
## Clean up
setwd("..")
