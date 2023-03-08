# Purpose:
#  Program for solving the system of coupled differential eqns
#  d_r^2 phi + 1/r dphi/dr - 1/r**2 (n + h a/r)**2 phi - m**2 phi - lambda phi**3  = 0

#  d_r^2 a - 1/r d_r a - h n phi**2 - h**2 a phi**2 = 0 
# with  boundary conditions
#    phi(0) = 0, phi(inf) = v, a(0) = 0, a(inf) = -n/h


import numpy as np
from scipy.integrate import solve_bvp
import scipy.special as spl
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

#mpl.rcParams.update({'font.size': 15})
# DEFINE
# Y = {phi, phi' ,a ,a' } = {y0,y1,y2,y3,y4,y5}
# Y'= {phi',phi'',a',a''} = {y1, -1/r y1 + ( n**2 + 2 h n a/r + h**2  a**2/r**2)y0/r**2 - mu**2  y0 + 2 lambda  y0**3 + kappa y0 y2**2,
#                                     y3, -1/r y3 + (n'**2 + 2 h'n'a/r + h'**2 a**2/r**2)y2/r**2 - mu'**2 y2 + 2 lambda' y2**3 + kappa y2 y0**2,
#                                     y5,  1/r y5 + h n y0**2 + h**2 y4 y0**2 + h' n' y2**2 + h'**2 y4 y2**2}
#

beta = 3.5
v    = 0.05
lmb  = beta/2#1.0
n    = 3.0
h   = 1#(2*lmb/beta)**0.5
a_inf = -n/h
msq  = -lmb * v**2 

def der(r,y):
    return np.vstack((y[1], -y[1]/r + y[0]*(n + h*(y[2]))**2/(r**2) +  msq*y[0] + lmb*(y[0])**3,\
                      y[3], y[3]/r  + h*n*(y[0])**2 + h**2*y[2]*(y[0])**2))

def bc(ya,yb):
    return np.array([ya[0],yb[0] - v,ya[2],yb[2] - a_inf])

# define mesh
r0 = 1e-3
rf = 100.0
N = 400
r = np.linspace(r0,rf,N)
y = np.zeros((4,r.size))

y_guess = np.zeros(r.size)
y2_guess = np.zeros(r.size)
dy_guess = np.zeros(r.size)
dy2_guess = np.zeros(r.size)

e = np.zeros(r.size)

datafile = open("local_profile.dat",'w')

for i in range(N):
    y_guess[i]  = v*(-math.exp(-r[i]**2)+1)
    y2_guess[i] = a_inf*(1+math.exp(-r[i]**2)) 
    dy_guess[i] = 2*v*r[i]*math.exp(-r[i]**2)
    dy2_guess[i]= 2*a_inf*r[i]*math.exp(-r[i]**2)#
    
    # initial guess
y[0,:] = y_guess   #  phi
y[1,:] = dy_guess  # dphi
y[2,:] = y2_guess  #  xi
y[3,:] = dy2_guess # dxi         
    #solve
res = solve_bvp(der,bc,r,y,tol = 0.001)
if res.status != 0:
  print("Warning: res.status is %d " % ( res.status))
    # plots 
y_plot   = res.sol(r)[0]
dy_plot  = res.sol(r)[1]
y2_plot  = res.sol(r)[2]
dy2_plot = res.sol(r)[3]

for i in range(N-1):
  e[i] = 0.5 * (dy_plot[i])**2 + 0.5 * n**2 * (1.0 + h * y2_plot[i])**2 * (y_plot[i])**2/(r[i])**2 \
  + 0.5 * msq * (y_plot[i])**2 + 0.25 * lmb * (y_plot[i])**4 \
  + 0.5 * (dy2_plot[i]/r[i])**2 - 0.5 * msq * v**2 - 0.25 * lmb * v**4
  print(r[i],y_plot[i],y2_plot[i],e[i], file = datafile)

fig = plt.figure(figsize=(8,5))    
ax1 = fig.add_axes([0.10,0.10,0.70,0.85])

ax1.plot(r,y_plot,r,y2_plot)
    
plt.text(90,0.7,'$\phi(r)$',fontsize = 16)
plt.text(90,-3.8,'$a(r)$',fontsize = 16)
#plt.text(8,0.05,'$\epsilon$',fontsize = 16)
plt.xlabel('$r$',fontsize =16)
plt.title("$n = $%d, $\\beta = $%3.1f, $\lambda = $%d,$v = $ %d"%(n,beta,lmb,v))
plt.show()
  
 
