`norm.gauss.hermite` <-
function(n){
# modify coefficients for standard normal
	rlist <- gauss.quad(n,kind="hermite")
	r <- cbind(rlist$nodes,rlist$weights)	
	r[,1] <- r[,1]*sqrt(2)
	r[,2] <- r[,2]/sqrt(pi)
	return(r)
}

