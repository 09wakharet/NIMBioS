#### Libraries
rm(list= ls())
library(graphics)
library(deSolve)
library(GA)

### Mathematical Model 
sir.single.outbreak = function(t,x,model.parameters){
  beta = model.parameters[1]
  gamma = model.parameters[2]
  N = x[1] + x[2] + x[3]
  eqns = vector(mode="numeric",length = length(x))
  eqns[1] = -beta*x[1]*x[2]/N
  eqns[2] = beta*x[1]*x[2]/N - gamma*x[2]
  eqns[3] = gamma*x[2]
  return(list(eqns))
}

### Initial conditions
x0.vec = rep(0,3)
x0.vec[1:2] = c(1e+3-1,1) 

### Time values
tspan = seq(0,6,by=5e-1)

### True parameter values
b0 = 5.5
g0 = 1.1
par0 = c(b0,g0)

### True data without noise
xt = lsoda(y=x0.vec,times=tspan, func=sir.single.outbreak, parms=par0)

#### Synthetic data with noise
tvals = tspan
xvals = xt[,3] + rnorm(length(xt[,3]),mean=0,sd=sqrt(100)) 





###############OLS WITH GLOBAL OPTIMIZATION

### Sum of squares
fit.ols.sir = function(theta){
  time.grid.fine = seq(tvals[1],tvals[length(tvals)],by=1e-2)
  out = lsoda(y=x0.vec,times=time.grid.fine,func=sir.single.outbreak,parms=theta);
  if(is.nan(out) || length(out[,2]) < length(time.grid.fine)){
    sum.squared.residuals = 1e+100
  }else{
    f.model = approx(time.grid.fine,out[,3],tvals)$y
    sum.squared.residuals = sum((f.model -xvals)^2)
  }
  return(sum.squared.residuals)
}

### Optimization: minimizing sum of square residuals with GA
lower.bounds = c(1,0.5)
upper.bounds = c(10,2.5)

fit.sir.ga = ga(type = "real-valued",
                fitness = function(q)-fit.ols.sir(q),
                maxiter = 5,
                min = lower.bounds,
                max = upper.bounds,
                popSize = 5e+3,
                parallel = TRUE)

parameter.estimates2 <- fit.sir.ga@solution

### plot dataset and bestfit solution

bestfit.solution = lsoda(y = x0.vec,times = seq(tvals[1],tvals[length(tvals)],by=1e-2),
                         func = sir.single.outbreak,parms = parameter.estimates2)

matplot(bestfit.solution[,1],bestfit.solution[,3],
        type="l",lwd=2,
        main="Prevalence Data and Best Fit Solution",ylab="I(t)",xlab="Time t")
points(tvals,xvals,col="dark green",pch=16,cex=2)


par0
parameter.estimates2
