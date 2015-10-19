import numpy as np
import astropy.constants as const 
import astropy.units as units
from scipy.integrate import quad
import cmath


def planck_function(nu, T):
  constant = 2.*h*nu**3 / c**2
  return constant / (np.exp(h*nu/k_B/T) - 1.)



def integrand(nu, T):
  return np.pi / h / nu * planck_function(nu, T) * alpha_nu* (nu_0/nu)**3






def integrand2(nu, T):
  
  return np.pi**2*4*r_0**2 * planck_function(nu,T)
  
  
def calculate_new_x(a,b,c):
  #solve ax**2+ bx + c = 0 and choose maximum
  
  #calculate discriminant
  d = (b**2) - (4*c)
  
  #calculate solutions
  sol1 = (-b-np.sqrt(d))/(2*a)
  sol2 = (-b+np.sqrt(d))/(2*a) 
  
  print sol1, sol2
  return max(sol1,sol2)
  
  
  
  
  
  









# user defined constants for this task

T_eff = 40e+3 # units.K
R = 20

alpha_nu = 6.3e-18 # units.cm**2
alpha = 2.59e-13 # units.cm**3 / units.s
n_H = 10 # units.cm**3
nu_0 = 3.2871e+15 # units.Hz

c = const.c.cgs.value 
h = const.h.cgs.value
k_B = const.k_B.cgs.value



#-----------------
# 
#-----------------


steps = 10000.
x_0 = 1.0
r_0 = R * 69550800000.0 # * units.cm
r_max = 100 * 3.085677581467192e+18

integral = quad(integrand2, 0. , nu_0,args=(T_eff))

#print integral

delta_r = (r_max-r_0)/steps

integral = quad(integrand, nu_0 , nu_0*1e+3,args=(T_eff))

A =  r_0**2 / alpha / n_H * integral[0]




x = np.zeros(int(steps+1))
r = np.zeros(int(steps+1))

#-----------------
# Iterative loop
#-----------------

#setting up the iteration
x[0] = x_0
r[0] = r_0

#begin iterative loop
for i in range(int(10)):
  
  iterative_sum = 0.
  
  for j in range(i+1):
    iterative_sum = (1.-x[j]) 
    
  
  f = 1./ (r[i]**2) * np.exp(-alpha_nu*n_H*delta_r*iterative_sum)
  #print x[i], f, r[i], A
  #print f
  

  
  x[i+1] = calculate_new_x(1.,A*f,A*f)
  
  #print x[i+1]
  #print x[i+1]
  r[i+1] = r[i] + delta_r
