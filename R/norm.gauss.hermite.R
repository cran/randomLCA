#*************************************************************************
#Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are
#met:
#
#- Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#- Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer listed
#  in this license in the documentation and/or other materials
#  provided with the distribution.
#
#- Neither the name of the copyright holders nor the names of its
#  contributors may be used to endorse or promote products derived from
#  this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES LOSS OF USE,
#DATA, OR PROFITS OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#*************************************************************************


#*************************************************************************
#Computation of nodes and weights for a Gauss-Hermite quadrature formula
#
#The  algorithm  calculates  the  nodes  and  weights  of the Gauss-Hermite
#quadrature  formula on domain  (-infinity, +infinity) with weight function
#W(x)=Exp(-x*x).
#
#Input parameters:
#    n   –   a required number of nodes.
#            1 <= n <= 190.
#
#Output parameters:
#    x   -   array of nodes.
#            Array whose index ranges from 0 to N-1.
#    w   -   array of weighting coefficients.
#            Array whose index ranges from 0 to N-1.
#
#The algorithm was designed by using information from the QUADRULE library.
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
        repeat
        {
            p2 <- 0
            p3 <- pipm4
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
