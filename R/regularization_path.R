#==========================================================#
# LASSO REGULARIZATION PATH solved with NETWORK FLOW
#==========================================================#

# *** Solve the L1 penalized poisson regression with network flow algorithm. 
# A dichotomy is performed in order to give a set of solution with a smoothly increasing nulber of active isoforms
# For each solution of the prenalized regression, a refit is performed *** 
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
      print(sprintf('the maximum value for lambda is %f', max_lambda_val))
   }
   if(max_lambda_val <= 0){
      max_lambda_val <- 1e-2
   }
   if(max_lambda_val > 1e10){
      max_lambda_val <- 1e10
   }
   #------ miscellaneous ------#
   param$lambda <- max_lambda_val
   dg.poisson.res <- duality_gap_poisson_graph(y=count, paths=matrix(0, 0, 0), 
                                               beta=matrix(0, 0, 0), graph=graph, param=param)
   primal0 <- dg.poisson.res$primal
   tmp0 <- max(count / param$loss_weights, na.rm=TRUE)
   # this field with be modified in-place, do not touch
   # it allows the flow algorithm to warm restart between two lambdas

   path.set <- list()
   beta.refit <- list()
   loss.set <- vector()
   size.set <- vector()

   if(verbose==1){ 
      print('REGULARIZATION PATH calculation ...')
   }

   ##================================
   ## Poisson: Penalized Regression
   ##================================
   # upper value
   lambda <- max_lambda_val
   max_capacity <- min(max(primal0/lambda,0), tmp0)
   epsilon_flow <- max(max_capacity*1e-9, 1e-10)
   if(fast_guess){
      paths <- candidate
      Beta <- betacandidate
      Z <- as.matrix(spams_flipflop.multLeftDiag(paths,param$loss_weights))
      Beta <- spams_flipflop.solverPoisson(y=count, X=as(Z,'dgCMatrix'), beta0=Beta, 
                                  weights=matrix(param$lambda, ncol(Z), 1),
                                  delta=param$delta, max_iter=iterpoisson, tol=tolpoisson)
      ind <- which(Beta != 0)
      if(length(ind) < ncol(paths) & ncol(paths)>1){ 
         paths <- paths[, ind]
         Beta <- Beta[ind]
      }
      paths <- as.matrix(paths)
      dg.res <- duality_gap_poisson_graph(y=count, paths=paths, beta=Beta, graph=graph, param=param)
      duality_gap <- dg.res$duality_gap
      primal <- dg.res$primal
      pc.res <- list()
      pc.res$path <- paths
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
      Beta <- matrix(apply(as.matrix(paths), 2, max), ncol=1)
      paths <- matrix(as.double(paths != 0),nrow=nrow(paths),ncol=ncol(paths))
      if(verbosepath){
         dg.res <- duality_gap_poisson_graph(y=count, paths=paths, beta=Beta, graph=graph, param=param)
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
   path.set[[1]] <- pc.res$path
   path.set[[1]][path.set[[1]]>0] <- 1
   beta0 <- matrix(apply(as.matrix(pc.res$path), 2, max), ncol=1)
   size.set[1] <- ncol(path.set[[1]])
   ## REFIT
   if(ncol(path.set[[1]])>0){
      Z <- as.matrix(spams_flipflop.multLeftDiag(paths,param$loss_weights))
      beta.refit[[1]] <- spams_flipflop.solverPoisson(y=count, X=as(Z,'dgCMatrix'), 
                                             beta0=beta0, weights=rep(0,length(beta0)), 
                                             delta=delta, max_iter=iterpoisson, tol=tolpoisson)
      loss.set[1] <- loss.ll(count, Z%*%beta.refit[[1]], noise.var=delta)
   }
   if(ncol(path.set[[1]])==0){
      beta.refit[[1]] <- 0
      loss.set[1] <- loss.ll(count, rep(0,n.nodes), delta)
   }
   lambdacase <- vector()
   sizecase <- vector()
   lambdacase[1] <- lambda
   numiso <- ncol(pc.res$path)
   sizecase[1] <- numiso
   lambda <- lambda/10

   while(lambda > 1e-3 && numiso < max_isoforms){
      # Compute isoforms given penalization lambda
      param$lambda <- lambda
      max_capacity <- min(max(primal0/lambda,0), tmp0)
      epsilon_flow <- max(max_capacity*1e-9, 1e-10)
      if(fast_guess){
         paths <- candidate
         Beta <- betacandidate
         Z <- as.matrix(spams_flipflop.multLeftDiag(paths,param$loss_weights))
         Beta <- spams_flipflop.solverPoisson(y=count, X=as(Z,'dgCMatrix'), 
                                     beta0=Beta, weights=matrix(param$lambda, ncol(Z), 1),
                                     delta=delta, max_iter=iterpoisson, tol=tolpoisson)
         ind <- which(Beta != 0)
         if(length(ind) < ncol(paths) & ncol(paths)>1){
            paths <- paths[, ind]
            Beta <- Beta[ind]
         }
         paths <- as.matrix(paths)
         dg.res <- duality_gap_poisson_graph(y=count, paths=paths, beta=Beta, graph=graph, param=param)
         duality_gap <- dg.res$duality_gap
         primal <- dg.res$primal
         pc.res <- list()
         pc.res$path <- paths
         if(verbosepath){
            cat(sprintf('first guess: num isoformes: %d, lambda: %f, primal: %f, 
                        duality-gap: %f\n', ncol(paths), param$lambda, primal, duality_gap))
         }
      }
      if(!fast_guess || duality_gap > 0.001) {
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
         Beta <- matrix(apply(as.matrix(paths), 2, max), ncol=1)
         paths <- matrix(as.double(paths != 0),nrow=nrow(paths),ncol=ncol(paths))
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
      sizemean <- ncol(pc.res$path)
      ind <- sizemean + 1 - min(sizecase, na.rm=TRUE)
      if(ind <=0){
         break
      }
      if(is.na(lambdacase[ind])==TRUE | (is.na(lambdacase[ind])==FALSE & lambdacase[ind]>lambda)){
         sizecase[ind] <- sizemean
         lambdacase[ind] <- lambda
         path.set[[ind]] <- pc.res$path
         path.set[[ind]][path.set[[ind]]>0] <- 1
         beta0 <- matrix(apply(as.matrix(pc.res$path), 2, max), ncol=1)
         ## REFIT
         if(ncol(path.set[[ind]])>0){
            Z <- as.matrix(spams_flipflop.multLeftDiag(paths, param$loss_weights))
            beta.refit[[ind]] <- spams_flipflop.solverPoisson(y=count, X=as(Z,'dgCMatrix'), 
                                                     beta0=beta0, weights=rep(0,length(beta0)), 
                                                     delta=delta, max_iter=iterpoisson, tol=tolpoisson)
            loss.set[ind] <- loss.ll(count, Z%*%beta.refit[[ind]], noise.var=delta)
            size.set[ind] <- length(which(beta.refit[[ind]]>0))
         }
         if(ncol(path.set[[ind]])==0){
            beta.refit[[ind]] <- 0
            loss.set[ind] <- loss.ll(count, rep(0,n.nodes), delta)
            size.set[ind] <- 0
         }
      }
      if(sum(is.na(sizecase)) > 0){	
         lambdacase[ind] <- lambda
         indna <- which(is.na(sizecase))[1]
         l1 <- lambdacase[indna-1]
         l2 <- max(lambdacase[indna:length(lambdacase)], na.rm=TRUE)
         if(abs(l2-l1) < 1e-5){
            break
         }
         lambdaold <- lambda
         lambda <- mean(c(l1,l2))
         if(lambdaold==lambda){
            break
         }
         numiso <- sizecase[indna-1] 
      } else{
         lambda <- min(lambdacase, na.rm=TRUE)/10
         numiso <- sizemean
      }
   } # End of regularization path

   na.final <- which(is.na(sizecase))
   if(sum(na.final)>0){
      beta.refit <- beta.refit[-na.final]
      path.set <- path.set[-na.final]
      size.set <- size.set[-na.final]
      loss.set <- loss.set[-na.final]
   }

   return(list(beta.refit=beta.refit, path.set=path.set, size.set=size.set, loss.set=loss.set))

}
