#
# logistic regression using our own IRWLS algorithm
# (c) Pedro Delicado
#
logistic.IRWLS <- function(x,y,weights.out=1,x.new=x,
                           max.iter=10,eps.beta=1e-5,
                           plts=FALSE){
  if (plts&(dim(as.matrix(x))[2]>1)){
    plts<-FALSE
    warning("Plots are valid only when dim(x)[2]==1")
  }
  # Step 0
  stop.rule <- FALSE
  iter <- 0
  beta.0 <- coef(lm(y~x)) 
  
  while (!stop.rule){
    iter <- iter + 1 
    # step 1
    p.lt <- p.from.beta(beta.0,x)
    p <- p.lt$p
    lt <- p.lt$lt
    ps.e <- (y-p)/(p*(1-p))
    z <- lt + ps.e 
    wt <- p*(1-p) *weights.out
    
    if (plts){
      op<-par(mfrow=c(1,2))
      plot(x,y,cex=8*wt)
      lines(x,p,col=2)
      plot(x,z,cex=8*wt)
      lines(x,lt,col=2)
      par(op)
    }
    
    lm.1 <- lm(z~x,weights = wt) 
    beta.1 <- coef(lm.1)
    
    # checking stop rules
    if ((iter>=max.iter)|(sum((beta.1-beta.0)^2)<eps.beta)){
      stop.rule<-TRUE
    } else {
      beta.0 <- beta.1
    }
  }
  
  aux <- summary(lm.1)
  p.lt <- p.from.beta(beta.1,x)
  p <- p.lt$p
  lt <- p.lt$lt
  se.beta <- diag(aux$cov.unscaled)^.5
  null.dev <- sum(-2*dbinom(y,1,mean(y),log=TRUE))
  resid.devi <- sum(-2*dbinom(y,1,p,log=TRUE))
  
  return(list(coefficients=beta.1, se.coef=se.beta,
              fitted.values=p,linear.predictors=lt,
              predicted.values=p.from.beta(beta.1,x.new)$p,
              residual.deviance=resid.devi,
              iter=iter)
         )
} 

# This is an auxiliary function for logistic.IRWLS
#
# p.from.beta(beta,x) 
# computes the vector of probabilities p, with elements
# p_i=\hat{Pr}(Y_i=1|X=x_i), i=1,...,n,
# where x_i is the ith row of matrix "x",
# a logistic link is assumed, and "beta" is the 
# estimated vector of coefficients.
p.from.beta <- function(beta,x){
  lin.term <- beta[1] + x %*% beta[-1]
  e.lt <- exp(lin.term)
  p <- e.lt/(1+e.lt)
  return(list(p=p,lt=lin.term))
}


######
#
# Description:
# Nonparametric logistic regression using spline smoothing
# fitted by a Penalized IRWLS algorithm.
# (c) Pedro Delicado
# 
# In this function there are essentially two changes 
# with respect to logistic.IRWLS:
# 
# 1. The linear model fit in line 
#         lm.1 <- lm(z~x,weights = wt) 
#    of logistic.IRWLS is replaced by a spline smoothing here:
#         spline.1 <- smooth.spline(x,z,w=wt,df=df,
#                         spar=spar, all.knots = all.knots, nknots = nknots) 
#
# 2. Calls to function p.from.beta(beta,x) in logistic.IRWLS are 
#    replaced here by calls to function p.from.theta(theta,x).
#
# Usage:
# logistic.IRWLS.splines(x=..., y=..., x.new=..., 
#              df=..., spar=..., all.knots = ..., 
#              nknots=..., max.iter=..., eps.beta=...,
#              plts=...)
# 
# Arguments:
# The amount of smoothing is controlled by the tuning parameter
# "df", that states the desired 'degrees of freedom' 
# (or 'equivalent number of parameters') for the fitted model.
# 
# All the parameters in this function (except "x.new" and "plts") have
# the same meaning that they have in function smooth.spline.
#
# The parameter "x.new" gives a new set of predicting variables 
# at which the estimated function must be evaluated, that is,
# the response variable will be predicted for each row of
# matrix "x.new" (using the estimated probability of Y being 1 
# given that X is equal to "x.new").
# 
# "plts" is logical value that indicates if plots have to be 
# done at each iteration of the algorithm (if plts==TRUE) 
# or not (if plts==FALSE).
#
# Value:
# An object with the following components:
# fitted.values, estimated probabilities p(x) of Y=1 given X=x 
#                for the training sample x.
# theta.x, estimation of theta(x)=log(p(x)/(1-p(x))) 
#                for the training sample x.
# df, equivalent number of parameters in the iterated calls to 
#                spline smoothing function spline.smoothing().
# predicted.values, estimated probabilities p(x.new) of Y=1 
#                given X=x.new for a new data set x.new.
# residual.deviance, residual deviance for the estimated model.
#
logistic.IRWLS.splines <- function(x,y,weights.out=1,x.new=x,
                           df=6,spar=NULL, 
                           all.knots = FALSE, nknots = .nknots.smspl,  
                           max.iter=10,eps.beta=1e-5,
                           plts=FALSE){
  if (plts&(dim(as.matrix(x))[2]>1)){
    plts<-FALSE
    warning("Plots are valid only when dim(x)[2]==1")
  }
  # Step 0
  stop.rule <- FALSE
  iter <- 0
  theta.0 <- fitted(lm(y~x)) 
  
  while (!stop.rule){
    iter <- iter + 1 
    # step 1
    p <- p.from.theta(theta.0)
    ps.e <- (y-p)/(p*(1-p))
    z <- theta.0 + ps.e 
    wt <- p*(1-p) *weights.out
    
    if (plts){
      op<-par(mfrow=c(1,2))
      plot(x,y,cex=8*wt)
      lines(x,p,col=2)
      plot(x,z,cex=8*wt)
      lines(x,theta.0,col=2)
      par(op)
    }
    
    spline.1 <- smooth.spline(x,z,w=wt,df=df,spar=spar,
                              all.knots = all.knots,
                              nknots = nknots) 
    theta.1 <- predict(spline.1,x=x)$y
    
    # checking stop rules
    if ((iter>=max.iter)|(sum((theta.1-theta.0)^2)<eps.beta)){
      stop.rule<-TRUE
    } else {
      theta.0 <- theta.1
    }
  }
  
  p <- p.from.theta(theta.1)
  resid.devi <- sum(-2*dbinom(y,1,p,log=TRUE))
  
  return(list(fitted.values=p,
              theta.x=theta.1,
              df=spline.1$df,
              predicted.values=p.from.theta(predict(spline.1,x=x.new)$y),
              residual.deviance=resid.devi)
  )
} 

# This is an auxiliary function for logistic.IRWLS.splines
#
# p.from.theta(theta,x) 
# computes the vector of probabilities p, with elements
# p_i=\hat{Pr}(Y_i=1|X=x_i), i=1,...,n,
# where x_i is the ith row of matrix "x",
# a logistic link is assumed, and "theta" is the vector 
# with the nonparametric estimations of 
# log( Pr(Y_i=1|X=x_i)/(1-Pr(Y_i=1|X=x_i)) ), i=1,...,n.
p.from.theta <- function(theta.x){
  p.x <- 1/(1+exp(-theta.x))
  return(p.x)
}
