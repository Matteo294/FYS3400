import numpy as np
from matplotlib import pyplot as plt 

eV_to_Joule = 1.602e-19 # Convert from electroVolt to Joule
xfact = 5.292e-11 # Length conversion factor for atomic units
Efact = 4.360e-18 # Energy conversion factor for atomic units
rho = 0.32e-10/xfact # Pauli repulsion parameter
lam = 1e3*eV_to_Joule/Efact # Lambda parameter

# Solve for f(x) = 0 -> Newton Raphson method
def NewtonRaphson(f, fdot, N, nsteps=1000, x=0):
    for _ in range(nsteps):
        x = x - f(x, N)/fdot(x, N)
    return x

# Partial Madelung constant
def alpha(N):
    s = 0
    for j in range(1, int(N/2)):
        s += (-1)**j/j
    return -2*s # Tested that gives 2ln(2) in the limit of large N

# I want to solve f(x) = 0 and this is f(x)
def f(x, N):
    return x**2*np.exp(-x/rho) - rho*alpha(N)/2/lam
# f'(x)
def fdot(x, N):
    return -2*alpha(N)/rho*np.exp(-x/rho) - rho*alpha(N)/2/lam

# Run from 0 to 100 particles
a = []
for n in range(0, 100):
    new_a = xfact*NewtonRaphson(f, fdot, n, x=1e-10/xfact) # solve and reconvert to normal units
    a.append(new_a)
    print(new_a)

plt.plot(range(100), a, '.')
plt.show()