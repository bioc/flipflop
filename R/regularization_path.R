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


regularization_path <- function(graph, count, param, max_isoforms, delta, fast_guess, iterpoisson, tolpoisson, verbosepath, verbose){

   n.nodes <- length(count)
   ## Compute the max lambda
   mlv.res <- max_lambda(y=count, graph=graph, 
                         delta=param$delta, loss_weights=param$loss_weights, 
                         pos=param$pos)
   max_lambda_val <- mlv.res$dual
   if(fast_guess){
      candidate <- mlv.res$path
      candidate <- matrix(as.double(candidate != 0), nrow=nrow(candidate), ncol=ncol(candidate))
      betacandidate <- matrix(0,1,1)
   }
   if(verbosepath){
      print(sprintf('the maximum value for lambda is %f', max_lambda_val))}
   if(max_lambda_val <= 0){ max_lambda_val <- 1e-2 }
   if(max_lambda_val > 1e10){ max_lambda_val <- 1e10 }
   
   #------ miscellaneous ------#
   param$lambda <- max_lambda_val
   dg.poisson.res <- duality_gap_poisson_graph(y=count, paths=matrix(0, 0, 0), 
                                               beta=matrix(0, 0, 0), graph=graph, param=param)
   primal0 <- dg.poisson.res$primal
   tmp0 <- max(count / param$loss_weights, na.rm=TRUE)
   # this field with be modified in-place, do not touch
   # it allows the flow algorithm to warm restart between two lambdas

   path.set <- list()
   beta.avantrefit <- list()
#   dualitygap.avantrefit <- vector()
#   loss.avantrefit <- vector()

   if(verbose==1){ 
      print('REGULARIZATION PATH calculation ...')}

   ##================================
   ## Poisson: Penalized Regression
   ##================================
   # upper value
   lambda <- max_lambda_val
   
   lambdacase <- vector()
   sizecase <- vector()
   lambdacase[1] <- lambda    # will be updated if it gives less shrinkage
   lambdasearch <- lambdacase # CHANGEMENT will be updated anyway

   numiso <- 1
   sizecase[1] <- ncol(candidate)
   path.set[[1]] <- candidate
   beta.avantrefit[[1]] <- betacandidate
#   dualitygap.avantrefit[1] <- dg.poisson.res$duality_gap
#   loss.avantrefit[1] <- primal0

   while(lambda > 1e-3 && numiso < max_isoforms){ 
      # Compute isoforms given penalization lambda
      param$lambda <- lambda # current lambda
      max_capacity <- min(max(primal0/lambda,0), tmp0)
      epsilon_flow <- max(max_capacity*1e-9, 1e-10)
      if(fast_guess){
         paths <- candidate
         Beta <- betacandidate ; 
         Z <- as.matrix(spams_flipflop.multLeftDiag(paths,param$loss_weights))
         Beta <- spams_flipflop.solverPoisson(y=count, X=as(Z,'dgCMatrix'), 
                                              beta0=Beta, weights=matrix(param$lambda, ncol(Z), 1),
                                              delta=delta, max_iter=iterpoisson, tol=tolpoisson)
         ind0 <- which(Beta != 0)
         if(length(ind0) < ncol(paths) & ncol(paths)>1){
            paths <- paths[, ind0, drop=F] ## CHANGEMENT
            Beta <- Beta[ind0,,drop=F] ## CHANGEMENT
         }
         dg.res <- duality_gap_poisson_graph(y=count, paths=paths, beta=Beta, graph=graph, param=param)
         duality_gap <- dg.res$duality_gap
         primal <- dg.res$primal
         if(verbosepath){
            cat(sprintf('first guess: num isoformes: %d, lambda: %f, primal: %f, 
                        duality-gap: %f\n', ncol(paths), param$lambda, primal, duality_gap))
         }
      }
      if(!fast_guess || duality_gap > 0.001){
         pc.res <- spams_flipflop.sepCostsPathCoding(count, 
                                                     graph, 
                                                     loss_weights=param$loss_weights,
                                                     max_capacity=max_capacity,
                                                     epsilon_flow=epsilon_flow,
                                                     lambda=lambda, 
                                                     delta=param$delta, 
                                                     tol=param$tol, 
                                                     regul=param$regul,
                                                     loss=param$loss,
                                                     pos=param$pos,
                                                     mode_decomposition=param$mode_decomposition)
         paths <- pc.res$path
         Beta <- matrix(apply(as.matrix(paths), 2, max), ncol=1)  ## OK attention à l'ordre !
         paths <- matrix(as.double(paths != 0),nrow=nrow(paths),ncol=ncol(paths)) ## OK
         if(verbosepath){
            dg.res <- duality_gap_poisson_graph(y=count, paths, Beta, graph, param)
            duality_gap <- dg.res$duality_gap
            primal <- dg.res$primal
            cat(sprintf('correction: num isoformes: %d, lambda: %f, primal: %f, 
                        duality-gap: %f\n', ncol(paths), param$lambda, primal, duality_gap))
         }
      }
      if(fast_guess && ncol(paths)>0){
         candidate <- paths
         betacandidate <- Beta
      }
      sizemean <- ncol(paths)
      ind <- sizemean + 1 - min(sizecase, na.rm=TRUE)
      if(ind <=0){ break }

      ##-----> case where you find a new set of solution or you have same solution size with a smaller lambda:
      if(is.na(lambdacase[ind])==TRUE | (is.na(lambdacase[ind])==FALSE & lambdacase[ind]>lambda)){
         sizecase[ind] <- sizemean
         lambdacase[ind] <- lambda # update the lambdas list
         lambdasearch[ind] <- lambda # update the search one
         path.set[[ind]] <- paths
         beta0 <- Beta 
         beta.avantrefit[[ind]] <- beta0
         #dualitygap.avantrefit[ind] <- duality_gap
         #loss.avantrefit[ind] <- primal
      }
      ##-----> case where you have same solution size with a bigger lambda:
      if((is.na(lambdacase[ind])==FALSE & lambdacase[ind]<lambda)){
         lambdasearch[ind] <- lambda # only update the search lambda
      }
      if(sum(is.na(sizecase)) > 0){ 	
         indna <- which(is.na(sizecase))[1] # position of first NA
         l1 <- lambdacase[indna-1]
         #l2 <- max(lambdacase[indna:length(lambdacase)], na.rm=TRUE)
         l21 <- max(lambdacase[indna:length(lambdacase)], na.rm=TRUE)
         l22 <- max(lambdasearch[indna:length(lambdasearch)], na.rm=TRUE)
         l2 <- max(c(l21, l22))
         if(abs((l2-l1)/l2) < 1e-3){  ## TODO remplacer break par oublier cette etape ... ?
            #print('ARREEEEEEEEET SUSPECT')
            break
         }
         lambdaold <- lambda
         lambda <- mean(c(l1,l2))
         if(lambdaold==lambda){
            #print('BREAAAAAAAAAAAAAAAAAK')
            break
         }
         numiso <- sizecase[indna-1]
      } else{
         lambdaold <- lambda
         lambda <- min(lambdacase, na.rm=TRUE)/10
         numiso <- max(sizecase) ## CHANGEMENT
         indre <- numiso + 1 - min(sizecase)
         candidate <- path.set[[indre]]  #  METTRE LE if fast_guess # CHANGEMENT
         betacandidate <- beta.avantrefit[[indre]] # CHANGEMENT
      }
   } # End of regularization path

   lambda.set <- lambdacase
   na.final <- which(is.na(sizecase))
   if(sum(na.final)>0){
      beta.avantrefit <- beta.avantrefit[-na.final]
      lambda.set <- lambdacase[-na.final]
      path.set <- path.set[-na.final]
      #loss.avantrefit <- loss.avantrefit[-na.final]
      #dualitygap.avantrefit <- dualitygap.avantrefit[-na.final]
   }

   ### REFIT 
   refit <- lapply(1:length(path.set), FUN=function(ind){
                   if(ncol(path.set[[ind]])>0){
                      Z <- as.matrix(spams_flipflop.multLeftDiag(path.set[[ind]], param$loss_weights))
                      beta.apos <- refitPoisson(y=count, X=as(Z,'dgCMatrix'),beta0= beta.avantrefit[[ind]], delta=delta,
                                                iterpoisson=iterpoisson, tolpoisson=tolpoisson, solver_refit='NEW', verbosepath=verbosepath)
                      loss.apos <- loss.ll(count, Z%*%beta.apos, delta)
                      size.apos <- length(which(beta.apos>0))
                      #dg.apos <- duality_gap_poisson(count, Z, beta.apos, delta)$dg
                   }
                   if(ncol(path.set[[ind]])==0){
                      beta.apos <- matrix(0,1,1)
                      loss.apos <- loss.ll(count, rep(0,n.nodes), noise.var=delta)
                      size.apos <- 0
                      #dg.apos <- duality_gap_poisson(count, matrix(0, n.nodes,1), 0, delta)$dg
                   }
                   if(verbosepath){
                      dg.apos <- duality_gap_poisson(count, Z, beta.apos, delta)$dg
                      cat(sprintf('refit: duality-gap: %f\n', dg.apos))
                   }
                   return(list(beta.apos=beta.apos, loss.apos=loss.apos, size.apos=size.apos)) })
   
   beta.refit  <- lapply(refit, FUN=function(k) k$beta.apos)
   loss.set  <- sapply(refit, FUN=function(k) k$loss.apos)
   size.set  <- sapply(refit, FUN=function(k) k$size.apos)
#   dualitygap.apresrefit <- sapply(refit, FUN=function(k) k$dg.apos)

   return(list(beta.refit=beta.refit, path.set=path.set, size.set=size.set, loss.set=loss.set, beta.avantrefit=beta.avantrefit, lambda.set=lambda.set
               ))

}
