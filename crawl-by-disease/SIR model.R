## Case Study for management paper
## Katie Scranton
## scranton.kt@gmail.com
## 9/27/2018

## SIR model

## We would like to find a region where the infected population cycles with long transient periods of almost no infected individuals
## in order to do this we need to likely add births and deaths (on the same time frame as waning immunity), then solve the 3 by 3 jacobian to find equilibrium values that give complex eigenvalues with real parts?
## would this be enough or do we need to add periodic forcing?

library(deSolve)

# Model parameters
## beta transmission
## gamma recovery
## eta vaccitation
theta.outbreak <- c(beta = 0.1,
           gamma = 0.01,
           eta = 0,
           nu=0)

theta.transient <- c(beta = 0.1,
                    gamma = 0.01,
                    eta = 0,
                    nu = 0.00001)

init.states <- c(S=0.9,
                 I=0.1,
                 R=0)

my.sir <- function(t,init.states,theta){
  with(as.list(c(init.states,theta)),{
    dS = -beta*S*I - eta*S + nu*R
    dI = beta*S*I - gamma*I
    dR = gamma*I + eta*S - nu*R
    
    list(c(dS,dI,dR))
  })
}

times <- seq(0,1000,by=0.1)

out.outbreak <- ode(y=init.states, times=times, func=my.sir, parms=theta.outbreak)
out.outbreak <- data.frame(out.outbreak)
out.transient <- ode(y=init.states, times=times, func=my.sir, parms=theta.transient)
out.transient <- data.frame(out.transient)

plot(out.transient$time,out.transient$I,type="l",ylim=c(0,1))
points(out.outbreak$time,out.outbreak$I,type="l",col="red")



