#==========================================================#
# LASSO REGULARIZATION PATH solved with NETWORK FLOW
# REFIT for no-shrinkage effect
#==========================================================#

# *** Solve the L1 penalized poisson regression with network flow algorithm. 
# A dichotomy is performed in order to give a set of solution with a smoothly increasing nulber of active isoforms
# For each solution of the penalized regression, a refit is performed *** 
# INPUT:   - splicing graph, vector of count
#          - maximum number of active isoforms in a solution (max_isoforms)
#          - solver parameters (all others arguments)
# OUTPUT:  - set of solution with for each the different paths (=structure of isoforms), associated abundances (beta), value of Poisson loss ans size (number of isoforms)


regularization_path_grouplasso_dich <- function(path.coll, beta.coll, count.samples, loss.weights, 
                                                delta, iterpoisson, tolpoisson, tolpoisson_refit, REFIT, verbosepath, verbose,
                                                n.samples, NN, len, max_isoforms, mc.cores){

   # compute the maximum lambda value
   Z <- as.matrix(spams_flipflop.multLeftDiag(path.coll, loss.weights))
   #beta.init <-  matrix(rep(beta.coll, n.samples), ncol=n.samples)
   beta.init <-  matrix(0, nrow=ncol(path.coll), ncol=n.samples) # start with 0 as you want 0 paths with max.lambda
   # max lambda val ? 
   grad1 <- 1.0 - count.samples/delta
   all.l2norm <- sapply( 1:ncol(Z), FUN=function(jj){
                        grad2 <- grad1*Z[,jj]
                        Tvec <- apply(grad2, 2, sum)
                        return(sqrt(sum(Tvec*Tvec))) })
   max.lambda <- max(all.l2norm)

   if(verbosepath) print(sprintf('the maximum value for lambda is %f', max.lambda))
   if(max.lambda <= 0){ max.lambda <- 1e-2 }
   #if(max.lambda > 1e10){ max.lambda <- 1e10 }

   path.set <- list()
   beta.avantrefit <- list()
   #   dualitygap.avantrefit <- vector()
   #   loss.avantrefit <- vector()

   if(verbose==1) print('REGULARIZATION PATH calculation ...')

   ##================================
   ## Poisson: Penalized Regression
   ##================================
   # upper value
   lambda <- max.lambda + 1000 # in order to be sure to start with 0 paths and not already 1.
   lambdacase <- vector()
   sizecase <- vector()
   lambdacase[1] <- lambda    # will be updated if it gives less shrinkage
   lambdasearch <- lambdacase # CHANGEMENT will be updated anyway
   beta.avantrefit <- list()

   beta.avantrefit[[1]] <- grouplassoPoisson(y=count.samples, X=as(Z,'dgCMatrix'), beta0=beta.init, delta=delta, 
                                             iterpoisson=iterpoisson, tolpoisson=tolpoisson, solver_refit='NEW', verbosepath=verbosepath, 
                                             n.samples=n.samples, lambda=lambda)

   numiso <- sum(apply(beta.avantrefit[[1]],1,sum)>0)
   sizecase[1] <- numiso
   #   dualitygap.avantrefit[1] <- dg.poisson.res$duality_gap
   #   loss.avantrefit[1] <- primal0

   lambda <- max.lambda
   beta.init <- beta.avantrefit[[1]]
   while(lambda > 1e-3 && numiso < max_isoforms){ # FAIRE UN FAST_GUESS ??  
      # compute isoforms given penalization lambda
      Beta <- grouplassoPoisson(y=count.samples, X=as(Z,'dgCMatrix'), beta0=beta.init, delta=delta, 
                                iterpoisson=iterpoisson, tolpoisson=tolpoisson, solver_refit='NEW', verbosepath=verbosepath, 
                                n.samples=n.samples, lambda=lambda)
      sizemean <- sum(apply(Beta,1,sum)>0)
      ind <- sizemean + 1 - min(sizecase, na.rm=TRUE)
      if(ind <=0){ break }

      #-----> case where you find a new set of solution or you have same solution size with a smaller lambda:
      if(is.na(lambdacase[ind])==TRUE | (is.na(lambdacase[ind])==FALSE & lambdacase[ind]>lambda)){
         sizecase[ind] <- sizemean
         lambdacase[ind] <- lambda # update the lambdas list
         lambdasearch[ind] <- lambda # update the search one
         beta.avantrefit[[ind]] <- Beta
         #dualitygap.avantrefit[ind] <- duality_gap
         #loss.avantrefit[ind] <- primal
      }

      #-----> case where you have same solution size with a bigger lambda:
      if((is.na(lambdacase[ind])==FALSE & lambdacase[ind]<lambda)){
         lambdasearch[ind] <- lambda # only update the search lambda
      }

      # decide how to choose the next lambda
      if(sum(is.na(sizecase)) > 0){	
         indna <- which(is.na(sizecase))[1] # position of first NA
         l1 <- lambdacase[indna-1]
         l21 <- max(lambdacase[indna:length(lambdacase)], na.rm=TRUE)
         l22 <- max(lambdasearch[indna:length(lambdasearch)], na.rm=TRUE)
         l2 <- max(c(l21, l22))
         lambdaold <- lambda
         lambda <- mean(c(l1,l2))
         if( abs((l2-l1)/l2) < 1e-3 | lambdaold==lambda ){  ## TODO remplacer break par oublier cette etape ... ?
            #print('BIZAAAAAAAAAAAAAAAAARE 1')
            sizecase[indna] <- sizecase[indna-1] # repeat artifical the previous solution if you cannot reach the current one
            lambdacase[indna] <- lambdacase[indna-1]
            lambdasearch[indna] <- lambdasearch[indna-1]
            beta.avantrefit[[indna]] <- beta.avantrefit[[indna-1]] # repeat artifical the previous solution if you cannot reach the current one
            if(sum(is.na(sizecase)) > 0){ 
               indna <- which(is.na(sizecase))[1] # position of first NA
               l1 <- lambdacase[indna-1]
               l21 <- max(lambdacase[indna:length(lambdacase)], na.rm=TRUE)
               l22 <- max(lambdasearch[indna:length(lambdasearch)], na.rm=TRUE)
               l2 <- max(c(l21, l22))
               if(lambdaold!=lambda){
                  lambda <- mean(c(l1,l2))
               } else {
                  lambda <- (3*l1+l2)/4
               }
            } else {
               lambda <- min(lambdacase, na.rm=TRUE)/10 # if there is no more NA after the replacement
            }
            #break
         }
         #if(lambdaold==lambda){
         #   print('BIZAAAAAAAAAAAAAAAAARE 2')
         #   break
         #}
         numiso <- sizecase[indna-1]
         beta.init <- beta.avantrefit[[indna-1]]
      } else {
         lambdaold <- lambda
         lambda <- min(lambdacase, na.rm=TRUE)/10
         numiso <- max(sizecase) ## CHANGEMENT
         #indre <- numiso + 1 - min(sizecase)
         #candidate <- path.set[[indre]]  #  METTRE LE if fast_guess # CHANGEMENT
         #betacandidate <- beta.avantrefit[[indre]] # CHANGEMENT
         beta.init <- beta.avantrefit[[length(beta.avantrefit)]]
      }
   } # End of regularization path

   lambda.set <- lambdacase
   na.final <- which(is.na(sizecase))
   if(sum(na.final)>0){
      beta.avantrefit <- beta.avantrefit[-na.final]
      lambda.set <- lambdacase[-na.final]
      sizecase <- sizecase[-na.final]
      #loss.avantrefit <- loss.avantrefit[-na.final]
      #dualitygap.avantrefit <- dualitygap.avantrefit[-na.final]
   }
   ind.keep <- which(!duplicated(sizecase))
   beta.avantrefit <- beta.avantrefit[ind.keep]
   lambda.set <- lambda.set[ind.keep]
   sizecase <- sizecase[ind.keep]

   ### REFIT 
   # refit the set of solutions and calculate global loss
   # !! refit only selected paths from the group lasso !!
   ind2refit <- lapply(1:length(beta.avantrefit), FUN=function(tt) which(apply(beta.avantrefit[[tt]],1,sum)>0))

   beta.refit <- replicate(length(beta.avantrefit), matrix(0, nrow=n.samples, ncol=ncol(path.coll)), simplify=FALSE)
   loss0 <- sum( apply( count.samples, 2, FUN=function(cc) loss.ll(cc, 0, delta) ) )
   loss.set <- rep(loss0, length(beta.avantrefit))
   size.set <- sizecase

   yhatall <- vector('list', length(beta.avantrefit))

   if(REFIT==TRUE){
      for(ii in 1:length(beta.avantrefit)){
         ind <- ind2refit[[ii]]
         yhatall[[ii]] <- matrix(0, nrow=nrow(count.samples), ncol=n.samples)
         if(length(ind)>0){
            path2fit <- path.coll[,ind,drop=F]
            beta.all <- beta.avantrefit[[ii]][ind,,drop=F]
            Z.samples <- lapply(1:n.samples, FUN=function(jj) as.matrix(spams_flipflop.multLeftDiag(path2fit, NN[jj]*len/1e9)))
            res.apos <- mclapply(1:n.samples, FUN=function(ind.c){
                               cc <- count.samples[,ind.c,drop=F]
                               Z <- Z.samples[[ind.c]]
                               # more appropriate initialization: either 0 vector or beta.avant refit (obtained with sum of counts)
                               # depending on l1 distance
                               if( sum(abs(cc - Z%*%beta.all[,ind.c,drop=F])) >= sum(abs(cc)) ){
                                  beta.init <- matrix(0,ncol(Z),1) # keep ????
                               } else {
                                  beta.init <- beta.all[,ind.c,drop=F]
                               }
                               beta.apos <- refitPoisson(y=cc, X=as(Z,'dgCMatrix'),beta0= beta.init, delta=delta,
                                                         iterpoisson=iterpoisson, tolpoisson=tolpoisson_refit, solver_refit='NEW', verbosepath=verbosepath)
                               loss.apos <- loss.ll(cc, Z%*%beta.apos, delta)
                               y.apos <- Z%*%beta.apos
                               return(list(beta.apos=beta.apos, loss.apos=loss.apos, y.apos=y.apos))}, mc.cores=mc.cores)

            beta.refit[[ii]][,ind] <- t( sapply( 1:n.samples, FUN=function(jj) return(res.apos[[jj]]$beta.apos) ) )
            loss.set[ii] <- sum( sapply( 1:n.samples, FUN=function(jj) return(res.apos[[jj]]$loss.apos) ) )
            size.set[ii] <- sum(apply(beta.refit[[ii]],2,sum)>0)
            yhatallii <- sapply( 1:n.samples, FUN=function(jj) return(res.apos[[jj]]$y.apos) )
            if(!is.matrix(yhatallii)) yhatallii <- matrix(yhatallii, ncol=n.samples)
            yhatall[[ii]] <- yhatallii
         }
      }
   } else {
      beta.refit <- lapply( beta.avantrefit, FUN=function(xx) t(xx))
      loss.set <- NULL
   }
   ## Calculate LOSS AS MERGE for testing
   if(REFIT==TRUE){
      loss.merge <- sapply(1:length(yhatall), FUN=function(jj) loss.ll( apply(count.samples,1,sum), apply(yhatall[[jj]],1,sum), delta ) )
   } else {
      loss.merge <- NULL
   }
   return(list(beta.refit=beta.refit, size.set=size.set, loss.set=loss.set, beta.avantrefit=beta.avantrefit, lambda.set=lambda.set, loss.merge=loss.merge
               )) 
}
