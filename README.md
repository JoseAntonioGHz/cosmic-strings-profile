# cosmic-strings-profile
This Python code solves numerically a system of non-linear coupled differential equations to obtain the profile of U(1)_Y' local cosmic strings.


# Purpose:
#  Program for solving the system of coupled differential eqns
#  d_r^2 phi + 1/r dphi/dr - 1/r**2 ( n**2 + 2 n h a/r + h**2  a**2/r**2) phi - m**2 phi - lambda phi**3 - kappa phi xi**2 = 0
#  d_r^2 xi  + 1/r dxi/dr  - 1/r**2 (n'**2 + 2 n'h'a/r + h'**2 a**2/r**2) xi  - m'**2 xi - lambda' xi**3 - kappa xi phi**2 = 0
#  d_r^2 a - 1/r d_r a - h**2 n phi**2 - h**2 a phi**2 - h' n' xi**2 - h'**2 a xi**2 = 0 
# with  boundary conditions
#    dphi(0) = 0, xi(0) = 0, phi(inf) = v, xi(inf) = v', a(0) = 0, a(inf) = -n/h = -n'/h'


phi = radial part of the higss boson
xi = radial part of the new higgs boson
a = auxiliary function of the U(1)_Y' gauge field 


It is believed that through the so-called Kibble mechanism topological defects could have formed in the early universe. We discuss a model
beyond the Standard Model that permits a type of topological defects called cosmic strings. 

In order to study cosmic string solutions, we first promote the global symmetry U(1)B−L to a local symmetry and add a new gauge coupling. 

The cancellation of gauge anomalies is achieved by adding a right-handed neutrino to each lepton generation. 

Moreover, a new Higgs field is added in order to give mass to the right-handed neutrino. 

Finally, we study the field equations of motion of the two Higgs fields and the gauge field, in
order to obtain the profiles of the cosmic strings.
