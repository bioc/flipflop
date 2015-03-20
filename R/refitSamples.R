refitSamples <- function(count.samples, path2fit, betasum, NN, len, delta, iterpoisson, tolpoisson, verbosepath, mc.cores){

   if(ncol(path2fit)>0) {
      Z.samples <- lapply(1:ncol(count.samples), FUN=function(jj) as.matrix(spams_flipflop.multLeftDiag(path2fit, NN[jj]*len/1e9)))
      
      beta.apos <- mclapply(1:ncol(count.samples), FUN=function(ind.c){
                          cc <- count.samples[,ind.c,drop=F]
                          Z <- Z.samples[[ind.c]]
                          # more appropriate initialization: either 0 vector or beta.avant refit (obtained with sum of counts)
                          # depending on l1 distance
                          if( sum(abs(cc - Z%*%betasum)) >= sum(abs(cc)) ){
                             beta.init <- matrix(0,ncol(Z),1) 
                          } else {
                             beta.init <- betasum
                          }
                          res <- refitPoisson(y=cc, X=as(Z,'dgCMatrix'),beta0= beta.init, delta=delta,
                                              iterpoisson=iterpoisson, tolpoisson=tolpoisson, solver_refit='NEW', verbosepath=verbosepath)
                          return(res)}, mc.cores=mc.cores)

      beta.apos <- t(do.call(cbind, beta.apos))

      if(verbosepath){
         dg.apos <- sapply(1:ncol(count.samples), FUN=function(ii) duality_gap_poisson(count.samples[,ii,drop=F], Z.samples[[ii]], beta.apos[ii,], delta)$dg)
         cat(sprintf('refit: duality-gap: %s\n', paste(dg.apos, collapse=' ')))
      }
   } else {
      beta.apos <- matrix(0,ncol(count.samples),1)
   }

   return(beta.refit=beta.apos)
}
