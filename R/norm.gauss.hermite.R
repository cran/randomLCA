#*************************************************************************
# HERMITE_COMPUTE computes a Gauss-Hermite quadrature rule.
#
#  Discussion:
#
#    The abscissas are the zeros of the N-th order Hermite polynomial.
#
#    The integration interval is ( -oo, +oo ).
#
#    The weight function is w(x) = exp ( - x*x );
#
#    The integral to approximate:
#
#      Integral ( -oo < X < +oo ) exp ( - x*x ) * F(X) dX
#
#    The quadrature rule:
#
#      Sum ( 1 <= I <= ORDER ) w(I) * F ( x(I) )
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#
#  Author:
#
#    FORTRAN77 original version by Arthur Stroud, Don Secrest
#    MATLAB version by John Burkardt
#    R conversion by ken Beath
#
#  Reference:
#
#    Arthur Stroud, Don Secrest,
#    Gaussian Quadrature Formulas,
#    Prentice Hall, 1966,
#    LC: QA299.4G3S7.
#
#  Parameters:
#
#    Input, n the order of the formula to be computed.
#
#    Output, r[,1], the Gauss-Hermite abscissas.
#
#    Output, r[,2], the Gauss-Hermite weights.
#
#*************************************************************************
`norm.gauss.hermite` <-
function(n){
    x <- rep(NA,n)    
    w <- rep(NA,n)
    
    pipm4 <- pi^(-0.25)
    for(i in 1:trunc((n+1)/2)) {
        if( i==1 )  r <-  sqrt(2*n+1)-1.85575*(2*n+1)^(-1/6)
        else  {
            if( i==2 ) 
                r <- r-1.14*(n^0.426)/r
             else
            {
                if( i==3 )  r <- 1.86*r+0.86*x[1]
                else{
                	if( i==4 )  r <- 1.91*r+0.91*x[2]
                    else  r <- 2*r+x[i-2]
                }
            }
        }
# inline rather than call hermite_root
		repeat
        {
            p2 <- 0
            p3 <- pipm4
# inline rather than call hermite_recur
            for(j in 1:n)
            {
                p1 <- p2
                p2 <- p3
                p3 <- p2*r*sqrt(2/j)-p1*sqrt((j-1)/j)
            }
            dp3 <- sqrt(2*(j-1))*p2
            r1 <- r
            r <- r-p3/dp3
       	 	if (!(abs(r-r1)>=.Machine$double.eps*(1+abs(r))*100)) break
        }
        x[i] <- -r
        w[i] <- 2/(dp3*dp3)
        x[n-i+1] <- -x[i]
        w[n-i+1] <- w[i]
     }
	r <- cbind(x,w)	
	r[,1] <- r[,1]*sqrt(2)
	r[,2] <- r[,2]/sqrt(pi)
	r[,2] <- r[,2]/sum(r[,2])
	return(r)
}
