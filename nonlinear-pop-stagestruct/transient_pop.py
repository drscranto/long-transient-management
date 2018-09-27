########
##
## Code for long transient model as published in:
## Morozov, Banerjee, Petrovskii, 2016 in Journal of Theoretical Biology
## Created on 10/28/2017
## Katherine Scranton scranton.kt@gmail.com
##
## all parameters have names following the paper:
## alpha, n, m, beta control the density dependent growth with Allee effects
## tau_A is the developmental delay
## we explicitly include delta juvenile mortality but we rename it d_J
## we rename adult mortality D as d_A
##
## we use the version of the model before it was dimensionless (eqn 11 in the reference)
## so that in the future we can explicitly change parameter values to match
## management interventions such as harvest of juveniles
##
########

from PyDDE import pydde
import numpy as np
from matplotlib.pyplot import *

## set global parameter values
tau_A = 45.
d_A = 0.5
d_J = 0.
beta = 1.
n_allee = 2
m_allee = 4
alpha = 1.

if d_A==0.45:
    init_A = 3.2
elif d_A==0.5:
    init_A = 1.9
else:
    init_A = 2.5

g = d_A*(beta**(1/m_allee))/(alpha*np.exp(-d_J*tau_A))
#print(g)

## define reproduction as a function of adult density
def rep(x):
    # is this beta term correct?
    return x*(alpha*x**(n_allee-1))/(1+(beta*x)**m_allee) 

## define the set of equations
def pop_grad(s,c,t):

    ## set a constant history: assume the pop was at a steady state at the init conditions
    lag_A = float(init_A)

    ## retrive the density of adults tau_A time units ago
    if t >= tau_A:
        lag_A = pydde.pastvalue(0,t-tau_A,0)

    ## calculate adult recruitment as:
    ## product of reproduction and survivorship over the delay
    R_A = rep(lag_A)*np.exp(-d_J*tau_A)
    
    dAdt = R_A - d_A*s[0] 
    
    return np.array([dAdt])

## solve from a set of initial conditions
init_cond = np.array([init_A])
dde_harvest = pydde.dde()
dde_harvest.dde(y=init_cond, times=np.arange(0.,4000.,1.), func=pop_grad, tol=1e-8,
           dt=0.1, nlag=1);

## plot
ion()
plot(dde_harvest.data[:,0], dde_harvest.data[:,1],label= str(g)[:6])
legend(loc='best')
xlabel('Time')
ylabel('Density')
savefig('pop_timeseries.pdf')
