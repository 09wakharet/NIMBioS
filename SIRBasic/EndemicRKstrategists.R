####################################################
#CUSTOM TWO HOST OPEN SIR MODEL (basis demography)
####################################################
##THESE CONDITIONS LEAD TO AN UNSTABLE SWITCH BETWEEN DYING OUT AND AN ENDEMIC STATE

SIRTwo <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    dS <- -bet*s.num*i.num - bet.g21*s.num*i2.num + alp *(s.num+i.num+r.num) - omeg*s.num
    dI <- bet*s.num*i.num-gam*i.num  + bet.g21*s.num*i2.num - omeg*i.num - omeg_star *i.num
    dR <- gam*i.num - omeg* r.num
    
    dS2 <- -bet.g2*s2.num*i2.num - bet.g12*s2.num*i.num + alp.g2 *(s2.num+i2.num+r.num) - omeg.g2 *s2.num
    dI2 <- bet.g2*s2.num*i2.num-gam.g2*i2.num + bet.g12*s2.num*i.num - omeg.g2*i2.num - omeg_star.g2 *i2.num
    dR2 <- gam.g2*i2.num - omeg.g2 * r2.num
    
    # Compartments and flows are part of the derivative vector
    # Other calculations to be output are outside the vector, but within the containing list
    list(c(dS, dI, dR, dS2, dI2, dR2))
  })
}

#beta transmission coefficients, gamma recovery, alpha birth, omega death
param <- param.dcm(bet =  c(1e-8,1e-7,1e-6,5e-6,10e-6), bet.g2 = 1.4e-6, bet.g21= 3e-6, bet.g12 = .2e-6, 
                   gam = .4, gam.g2 = .7,
                   omeg = .1, omeg.g2 = .01, 
                   omeg_star = .01, omeg_star.g2 = .01,
                   alp = .11, alp.g2 = .02) #balance av. contact rate between g1 and g2)?
init <- init.dcm(s.num = 1e6, i.num = 1, r.num = 0,
                 s2.num = 1e3, i2.num = 1, r2.num =0)
control <- control.dcm(nsteps = 150, dt = 0.02, new.mod = SIRTwo)

mod <- dcm(param, init, control)
head(as.data.frame(mod)[,1:3])

par(mfrow = c(1, 1))
#plot(mod, y = c("s.num","s2.num","i.num","i2.num","r.num","r2.num"), main = "Two Host Dynamics", legend = "full")
plot(mod, y = "i.num", main = "beta 1 senstivity analysis", legend = "full")

