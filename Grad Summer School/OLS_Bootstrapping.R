rm(list = ls()); ## clear memory
library(graphics) ## package for plots
library(deSolve)
library(GA) ## package for genetic algorithm optimization 

##### Mathematical model
logistic.model = function(t,x,paras){
  ## model parameters
  r = paras[1];
  K = paras[2];
  ## model equation
  xdot = r*x*(1-(x/K));
  return(list(xdot))
}
#### True parameter values
r0 = 5e-1;
K0 = 7e+2;
true.param.vector = c(r0,K0)
### True data without noise
timevalues0 = seq(0,10,by=0.5)
IC.true = 30;
xtrue = lsoda(y=IC.true, times=timevalues0,func=logistic.model,parms=true.param.vector)
### Synthetic data
tvals = timevalues0
xvals = xtrue[,2] + rnorm(length(timevalues0),mean=0,sd=25) 
### Sum of squared residuals
fit.ols = function(p){
  time.grid.fine = seq(tvals[1],tvals[length(tvals)],by=1e-2)
  out = lsoda(y=xvals[1],times=time.grid.fine,func=logistic.model,parms=p);
  if(is.nan(out) || length(out[,2]) < length(time.grid.fine)){
    sum.squared.residuals = 1e+100
  }else{
    f.model = approx(time.grid.fine,out[,2],tvals)$y
    sum.squared.residuals = sum((f.model -xvals)^2)
  }
  return(sum.squared.residuals)
}

### Optimization: minimizing sum of square residuals with GA

LowerBounds = c(0,300)
UpperBounds = c(1,1000)

fit.ga1 = ga(type="real-valued",
             fitness=function(x)-fit.ols(x),
             maxiter=2.5e+1, #increase this argument
             min=LowerBounds,
             max = UpperBounds,
             popSize=1e+3)

parameter.estimates1 <- fit.ga1@solution



### plot dataset and bestfit solution
numsolBestFit = lsoda(y=xvals[1],times=seq(tvals[1],tvals[length(tvals)],by=1e-2),
                      func=logistic.model,parms=parameter.estimates1)
matplot(numsolBestFit[,1],numsolBestFit[,2],type="l",lwd=2,
        main="Data and Best Fit Solution",ylab="x(t)",xlab="Time t")
points(tvals,xvals,col="red",pch=16,cex=2)


### Displaying carrying capacity
numsolBestFit2 = lsoda(y=xvals[1],times=seq(tvals[1],tvals[length(tvals)]+10,by=1e-2),
                       func=logistic.model,parms=parameter.estimates1)
matplot(numsolBestFit2[,1],numsolBestFit2[,2],type="l",lwd=2,
        main="Data and Best Fit Solution",ylab="x(t)",xlab="Time t")
points(tvals,xvals,col="red",pch=16,cex=2)


###parameter values
true.param.vector
parameter.estimates1
