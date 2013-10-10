##========================================
## FUNCTIONS FOR THE NETWORK FLOW SOLVER
##=========================================
spams_flipflop.multLeftDiag <- function(X,Y) {
#  XAt = matrix(rep(0,m * n),nrow = m,ncol = n)
  m=nrow(X)
  n=ncol(X)
  XY= matrix(c(0),nrow = m,ncol = n)
  multLeftDiag(X,Y,XY)
  return(XY)
}

spams_flipflop.fistaFlat <- function(Y,X,W0,return_optim_info = FALSE,numThreads =-1,max_it =1000,L0=1.0,
              fixed_step=FALSE,gamma=1.5,lambda1=1.0,delta=1.0,lambda2=0.,lambda3=0.,
              a=1.0,b=0.,c=1.0,tol=0.000001,it0=100,max_iter_backtracking=1000,
              compute_gram=FALSE,lin_admm=FALSE,admm=FALSE,intercept=FALSE,
              resetflow=FALSE,regul="",loss="",verbose=FALSE,pos=FALSE,clever=FALSE,
              log=FALSE,ista=FALSE,subgrad=FALSE,logName="",is_inner_weights=FALSE,
              inner_weights=c(0.),size_group=1,sqrt_step=TRUE,transpose=FALSE) {

  m = nrow(W0)
  n = ncol(W0)
#  W = matrix(rep(0,m * n),nrow = m,ncol = n)
  W = matrix(c(0),nrow = m,ncol = n)
#  optim_info = do.call(solver,c(list(Y,X,W0,W),params))
##  optim_info = .mycall('fistaFlat',c('Y','X','W0','W',params))
  optim_info = fistaFlat(Y,X,W0,W,numThreads ,max_it ,L0,fixed_step,gamma,lambda1,delta,lambda2,lambda3,a,b,c,tol,it0,max_iter_backtracking,compute_gram,lin_admm,admm,intercept,resetflow,regul,loss,verbose,pos,clever,log,ista,subgrad,logName,is_inner_weights,inner_weights,size_group,sqrt_step,transpose)
  if(return_optim_info == TRUE)
    return(list(W,optim_info))
  else
    return (W)
}

spams_flipflop.solverPoisson <- function(y,X,beta0,weights,delta, max_iter=500, tol=1e-4) {
   if(class(X) != 'dgCMatrix'){
      stop("X should be a sparse matrix")
   }
   beta=matrix(c(0),nrow= ncol(X),1)
   solverPoisson(y,X,beta0,beta,weights,delta=delta,max_iter,tol)
   return(beta)
}

spams_flipflop.solverPoissonFull <- function(y,X,beta0,weights,delta, max_iter=500, tol=1e-4) {
   beta=matrix(c(0),nrow= ncol(X),1)
   solverPoissonFull(y,X,beta0,beta,weights,delta=delta,max_iter,tol)
   return(beta)
}


spams_flipflop.sepCostsPathCoding <- function(alpha0, DAG, loss_weights, max_capacity=1e10, epsilon_flow=1e-10, prices=NULL,
                                     numThreads =-1, lambda=1.0, 
                                     regul="", loss="", pos=FALSE, tol=1e-5, delta=1e-3, mode_decomposition=1){

  if(class(prices) == 'dgCMatrix'){
    stop("sepCostsPathCoding : prices should not be a sparse matrix")
  }
  
  if (length(DAG) != 3) {
   stop("sepCostsPathCoding : DAG should be a list of 3 elements")
  }
  
  if(class(DAG[['weights']]) != 'dgCMatrix'){
    stop("sepCostsPathCoding : DAG[['weights']] should be a sparse matrix of class dgCMatrix")
  }
  
  start_weights = DAG[['start_weights']]
  stop_weights =  DAG[['stop_weights']]
  ir = DAG[['weights']]@i
  jc = DAG[['weights']]@p
  weights = DAG[['weights']]@x

  if(is.null(prices)){
    n <- length(start_weights)
    prices <- matrix(c(0), nrow=2*n+2, ncol=1)
  }
  
  alpha = matrix(c(0),nrow = nrow(alpha0),ncol = ncol(alpha0))
  path = sepCostsPathCoding(alpha0,alpha,weights, ir, jc, start_weights, stop_weights, max_capacity, epsilon_flow, prices,
    numThreads,lambda,tol,delta,loss_weights,regul,loss,pos,mode_decomposition)

  indptr = path[[1]]
  indices = path[[2]]
  data = path[[3]]
  shape = path[[4]]
  path.sp.mat = sparseMatrix(i = indices, p = indptr, x = data, dims = shape, index1 = FALSE)    

  return(list('alpha'=alpha, 'path'=path.sp.mat))
}

spams_flipflop.evalPathCoding <- function(alpha, DAG, numThreads =-1, lambda1=1.0, lambda2=0.,
                 intercept=FALSE, resetflow=FALSE, regul="",verbose=FALSE, precision=1e9,
                 pos=FALSE, clever=TRUE, eval=FALSE, eval_dual=FALSE, size_group=1, transpose=FALSE) {
  
  if (length(DAG) != 3) {
   stop("evalPathCoding : DAG should be a list of 3 elements")
  }
  
  if(class(DAG[['weights']]) != 'dgCMatrix'){
    stop("evalPathCoding : DAG[['weights']] should be a sparse matrix of class dgCMatrix")
  }
  
  start_weights = DAG[['start_weights']]
  stop_weights =  DAG[['stop_weights']]
  ir = DAG[['weights']]@i
  jc = DAG[['weights']]@p
  weights = DAG[['weights']]@x

  val <- c(0)
  ##val <- matrix(c(0),nrow = 1,ncol = 1)
  
  path = evalPathCoding(alpha, val,
    precision, weights, ir, jc, start_weights, stop_weights,
    numThreads,lambda1,lambda2,intercept,resetflow,regul,verbose,pos,clever,eval,eval_dual,size_group,transpose)

  indptr = path[[1]]
  indices = path[[2]]
  data = path[[3]]
  shape = path[[4]]
  path.sp.mat = sparseMatrix(i = indices, p = indptr, x = data, dims = shape, index1 = FALSE)    

  ##return(list('alpha'=alpha, 'path'=path.sp.mat, 'dual'=val))
  
  if(eval_dual){
    return(list('dual'=val, 'alpha'=alpha, 'path'=path.sp.mat))}
  else{
    return(list('alpha'=alpha, 'path'=path.sp.mat))}
}
