require(EpiModel)
require(deSolve)

SIRTwo <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    dS <- -bet*s.num*i.num - bet.g21*s.num*i2.num + alp * (1+sin(2*3.14159*t))*(s.num+i.num+r.num) - omeg*s.num
    dI <- bet*s.num*i.num-gam*i.num  + bet.g21*s.num*i2.num - omeg*i.num - omeg_star *i.num
    dR <- gam*i.num - omeg* r.num
    
    dS2 <- -bet.g2*s2.num*i2.num - bet.g12*s2.num*i.num + alp.g2 *(s2.num+i2.num+r.num) - omeg.g2 *s2.num
    dI2 <- bet.g2*s2.num*i2.num-gam.g2*i2.num + bet.g12*s2.num*i.num - omeg.g2*i2.num - omeg_star.g2 *i2.num
    dR2 <- gam.g2*i2.num - omeg.g2 * r2.num
    
    #anything with a g2 subscript referes to population 2    
    list(c(dS, dI, dR, dS2, dI2, dR2))
    #list(c(dS, dI, dR, dS2, dI2, dR2), prop.infected = i.num/(s.num+i.num+r.num))
  })
}

#beta transmission coefficients, gamma recovery, alpha birth, omega death
#at equilibrium in the absence of pathogen
dem_scale_factor = .6 #multiples all demographic processes
param <- param.dcm(bet = 3e-6, bet.g2 = 3e-6, 
                   bet.g21= 1e-7, bet.g12 = 1e-7, 
                   gam = 1e-6, gam.g2 = 1e-6,
                   omeg = dem_scale_factor, omeg.g2 = dem_scale_factor, 
                   omeg_star = dem_scale_factor*1e-2, omeg_star.g2 = dem_scale_factor*1e-2,
                   alp = dem_scale_factor, alp.g2 = dem_scale_factor) #balance av. contact rate between g1 and g2)?
init <- init.dcm(s.num = 1e6, i.num = 1, r.num = 0,
                 s2.num = 1e6, i2.num = 1, r2.num =0)
control <- control.dcm(nsteps = 10, dt = 0.02, new.mod = SIRTwo)

mod <- dcm(param, init, control)
head(as.data.frame(mod)[,1:3])

par(mfrow = c(1, 1))
#plot(mod, y = c("s.num","s2.num","i.num","i2.num","r.num","r2.num"), main = "Two Host Dynamics", legend = "full")
#plot(mod, y = c("s.num","i.num","r.num","prop.infected"), main = "beta 1 senstivity analysis", legend = "full")
plot(mod, y = c("s.num","s2.num"), main = "beta 1 senstivity analysis", legend = "full")



#######################
#TODO: mess with alpha and omega
#tip: if you add extra lines to the equation output outside the c() 
#(like the commented stuff), it'll give you extra columns you can plot
######################