##===================================================
## POISSON loss
##===================================================
loss.ll <- function(x, w, noise.var) 
{
	n.nodes <- length(x)
	sum(w - x*log(w + matrix(noise.var,nrow=n.nodes,ncol=1))) 
}
##===================================================

