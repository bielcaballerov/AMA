# locpolreg.R Local polynomial regression for estimating the 
#             regression function or its r-th derivative
#            
# Input: 
#      x,y  Observed data (two (n,1) vectors)
#      h    Smoothing parameter 
#      q    degree of the local polynomial to be fitted (default: 1)
#      r    order of the derivative to be estimate (Default: 0, the function)
#      tg   grid of values t where the estimated regression function 
#           is evaluated (default: x)
#      type.kernel "normal"  (Gaussian, default), 
#                  "epan"    (Epanechnikov) or 
#                  "rs.epan" (re-scaled Epanechnikov)
#                  "unif"    (Uniform Kernel in [-1,1])
#
# Output:  An object with two elements, 
#      mtg  Estimated values of the r-th derivative of the regression function at points in vector tg
#      S    The smoothing matrix
#
# Use:
#      result <- locpolreg(x,y,h,q,r,tg,type.kernel)
# 
locpolreg <- function(x,y,h=quantile(dist(x),.25),q=1,r=0,tg=NULL,type.kernel="normal",
                      nosubplot=FALSE,doing.plot=TRUE, ...){
   if (is.null(tg)){tg<-x}                  
   aux <- sort(tg,index.return=T)
   sorted.tg <- tg[aux$ix]
   sorted.tg.ix <- aux$ix

   n <- length(x);
   m <- length(tg);
   mtgr <- numeric(m);
   S <- matrix(0,nrow=m,ncol=n)

   for (i in seq(1,m)){
      aux <- kernel((x-tg[i])/h,type=type.kernel);
      Ih <- (aux>0);
      ni <- sum(Ih);     
      xh <- x[Ih]-tg[i];
      Dq <- matrix(1,nrow=ni,ncol=q+1);
      if (q>0){for (j in 1:q) Dq[,j+1] <- xh^j}
      Wx <- kernel(xh/h,type=type.kernel)/h;
      Wm <- Wx%*%ones(1,q+1);
      Dqq <- Wm*Dq;
      Si <- solve(t(Dq)%*%Dqq)%*%t(Dqq);
      beta <- Si%*%y[Ih];
      mtgr[i] <- factorial(r)*beta[r+1];
      S[i,Ih] <- Si[r+1,]
   }
  
   if (doing.plot){
      if (r==0){
        if (nosubplot) par(mfrow=c(1,1))
        plot(x,y,col="grey",...)
        lines(sorted.tg,mtgr[sorted.tg.ix],col=1,lwd=2)
      } 
      else{
         par(mfrow=c(2,1))
         aux <- locpolreg(x,y,h,q,0,tg,nosubplot=F,type.kernel,...)
         plot(sorted.tg,mtgr[sorted.tg.ix],type="n", 
              xlab="x",ylab="Estimated derivative")
         abline(h=0,col=4)
         lines(sorted.tg,mtgr[sorted.tg.ix],col=1,lwd=2)
      }
   }
return(list(mtgr=mtgr,S=S))
}

epan <- function(x){pmax(.75*(x+1)*(1-x))}
kernel <- function(x,type=c("normal","epan","rs.epan","unif")){
   switch(type[1],
          epan = pmax(.75*(x+1)*(1-x),0),
          rs.epan = pmax(.75*(x/sqrt(5)+1)*(1-x/sqrt(5))/sqrt(5),0),
          unif = as.numeric( (abs(x)<=1) )/2,
          dnorm(x))
}
ones <- function(n,m){matrix(1,nrow=n,ncol=m)}
