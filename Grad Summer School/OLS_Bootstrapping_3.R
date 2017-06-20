rm(list = ls()); ## clear memory
library(deSolve)
library(graphics)
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
timevalues0 = seq(0,20,by=0.5)
IC.true = 30;
xtrue = lsoda(y=IC.true, times=timevalues0,func=logistic.model,parms=true.param.vector)
### Synthetic data
tvals = timevalues0
xvals = xtrue[,2] + rnorm(length(timevalues0),mean=0,sd=100) 
xvals1 = xvals
matplot(tvals,xvals,col="blue",pch=16,cex=2, xlab="Time t",ylab="x(t)")

### Sum of squared residuals
fit.ols = function(p){
  time.grid.fine = seq(tvals[1],tvals[length(tvals)],by=1e-2)
  out = lsoda(y=IC.true,times=time.grid.fine,func=logistic.model,parms=p);
  if(is.nan(out) || length(out[,2]) < length(time.grid.fine)){
    sum.squared.residuals = 1e+100
  }else{
    f.model = approx(time.grid.fine,out[,2],tvals)$y
    sum.squared.residuals = sum((f.model -xvals)^2)
  }
  return(sum.squared.residuals)
}

### Optimization: minimizing sum of square residuals
initial.guess = c(0.25,350) ## initial guess
LowerBounds = c(0,300)
UpperBounds = c(1,1000)
fit.optim1 = optim(par=initial.guess,fn=fit.ols,method="L-BFGS-B",
                   lower=LowerBounds,upper=UpperBounds)
parameter.estimates1 <- fit.optim1$par

numsolBestFit = lsoda(y=IC.true,times=tvals,
                      func=logistic.model,parms=parameter.estimates1)
scaling.fac =  sqrt(length(xvals)/(length(xvals)-length(parameter.estimates1)))
std.res = scaling.fac*(xvals - numsolBestFit[,2])


number.bootstrap.samples = 20
boostrap.matrix = matrix(nrow = number.bootstrap.samples,ncol = length(parameter.estimates1))
boostrap.samplepoints.matrix = matrix(nrow = length(tvals),ncol = number.bootstrap.samples)
for(m in seq(1,number.bootstrap.samples,by=1)){
  xvals = numsolBestFit[,2] + sample(std.res,replace = TRUE) ## over-write xvals with bootstrap sample points
  fit.logistic.bootstrap = optim(par=initial.guess,fn=fit.ols,method="L-BFGS-B",
                                 lower=LowerBounds,upper=UpperBounds)
  boostrap.matrix[m,] = fit.logistic.bootstrap$par
  bestfit.sol.bootstrap = lsoda(y=IC.true,times = tvals,
                                func = logistic.model,
                                parms = fit.logistic.bootstrap$par)
  boostrap.samplepoints.matrix[,m] = bestfit.sol.bootstrap[,2]
}

matplot(tvals,boostrap.samplepoints.matrix,type="l",col="gray",lty=1,
        xlab="Time t",ylab="x(t)")

true.param.vector
estimates.bootstrap = parameter.estimates1
for(k in 1:dim(boostrap.matrix)[2]){
  estimates.bootstrap[k] = mean(boostrap.matrix[,k])
}
estimates.bootstrap
c.matrix = cov(boostrap.matrix)
sqrt(diag(c.matrix))
points(tvals,xvals1,col="blue",pch=16,cex=2)