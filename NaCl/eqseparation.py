import numpy as np
from matplotlib import pyplot as plt 
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

eV_to_Joule = 1.602e-19 # Convert from electroVolt to Joule
xfact = 5.292e-11 # Length conversion factor for atomic units
Efact = 4.360e-18 # Energy conversion factor for atomic units
rho = 0.32e-10/xfact # Rho parameter in the Pauli interaction lam*exp(-r/rho)
lam = 1e3*eV_to_Joule/Efact # Lambda parameter in the Pauli interaction lamba*exp(-r/rho)

# Solve for f(x) = 0 with the Newton Raphson method (in this case f(x) is the derivarive of the total energy U)
def NewtonRaphson(f, fdot, N, nsteps=1000, x=0):
    for _ in range(nsteps):
        x = x - f(x, N)/fdot(x, N)
    return x

# Partial Madelung constant
def alpha(N):
    s = 0
    for j in range(1, int(N/2)):
        s += (-1)**j/j
    return -2*s

# I want to solve f(x) = 0 and this is f(x)
def f(x, N):
    return x**2*np.exp(-x/rho) - rho*alpha(N)/2/lam # U'(x) = 0
# f'(x) used in the Newton-Raphson method
def fdot(x, N):
    return -2*alpha(N)/rho*np.exp(-x/rho) - rho*alpha(N)/2/lam

# Run from 0 to 100 particles
Nparticles = np.arange(0, 101, 2)
a = []
for n in Nparticles:
    new_a = xfact*NewtonRaphson(f, fdot, n, x=1e-10/xfact) # solve and reconvert to normal units
    a.append(new_a/1e-10)
    print("N atoms:", n, "Lattice constant:", new_a, "angstrom")

plt.plot(Nparticles, a, '.-', color='cornflowerblue', markerfacecolor='darkorange', markeredgecolor='black', markersize=10, linewidth=2.2)
plt.grid()
plt.xlabel('Number of atoms', fontsize=16)
plt.ylabel(r'Lattice constant [$\AA$]', fontsize=16)
plt.show()