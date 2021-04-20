import numpy as np 
from matplotlib import pyplot as plt
from scipy.optimize import fsolve 
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Relation between n=2N and ef, setting V=1 m=1 hbar=1 ef=1
Nelect1 = 1/3/np.pi**2 * np.power(2, 3/2)

# Numerical integral (finite sum) of the average number of fermions
def Nelect2(mu, T):
    x = np.linspace(0, 1000, 10000)
    dx = x[1] - x[0]
    return sum(DOS_fFD(x, mu, T) * dx)

# Integrand function
def DOS_fFD(eps, mu, T):
    return 1/2/np.pi**2 * 2**(3/2) * np.power(eps, 1/2) * fFD(eps, mu, T)

# Fermi-Dirac distribution function with kB=1
def fFD(eps, mu, T):
    return 1/(np.exp((eps-mu)/T) + 1)

# Equation f(x) = 0 to solve numerically
def eq(mu, T):
    return Nelect2(mu, T) - Nelect1

# Solver
T_vals = np.linspace(0.01, 10, 100) # from T=0.01Tf to T=10Tf
mu_vals = []
T_vals_zoom = np.linspace(0.01, 0.09, 100)
mu_vals_zoom = []
for T in T_vals:
    mu_vals.append(fsolve(eq, x0=1, args=(T)))
for T in T_vals_zoom:
    mu_vals_zoom.append(fsolve(eq, x0=1, args=(T)))

# Plot of mu(T)
plt.plot(T_vals, mu_vals, linewidth=2.0, color='black')
plt.xlabel(r'$T/T_F$', fontsize=18)
plt.ylabel('$\mu(T)$', fontsize=18)
plt.grid()
# Zoom in
fig = plt.figure(1)
ax_new = fig.add_axes([0.6, 0.6, 0.2, 0.2])
plt.plot(T_vals_zoom, mu_vals_zoom, linewidth=1.6, color='black')
plt.show()


''' Study of the Fermi-Dirac distribution '''
T_fFD = [0.01, 0.10, 0.50, 1.00, 1.20] # T vals for the Fermi-Dirac distribution study
#plot_colors = ['lightskyblue', 'turquoise', 'springgreen', 'red', 'blue']
plot_colors = ['blue', 'springgreen', 'turquoise', 'red', 'gold']
# Study at constant mu=mu(T=0)
#---------------------------------------------
def fFD_mu0(x, T, mu=1,):
    return 1/(np.exp((x-mu)/T) + 1)

eps_vals = np.linspace(0, 3, 10000)
for i in range(len(T_fFD)):
    FD_T = fFD_mu0(eps_vals, T_fFD[i])
    plt.plot(eps_vals, FD_T, color=plot_colors[i], label=r'T=%.2f $T_F \quad \mu=\epsilon_f$'%(T_fFD[i]))
#--------------------------------------------

# Study at varying mu=mu(T)
#-------------------------------------------
# mu vals for the Fermi-Dirac distribution study
mu_fFD = []
for T in T_fFD:
    mu_fFD.append(fsolve(eq, x0=1, args=(T)))

for i in range(len(mu_fFD)):
    plt.plot(eps_vals, fFD(eps_vals, mu_fFD[i], T_fFD[i]), '--', alpha=0.3, color=plot_colors[i], label=r'T=%.2f $T_F \quad \mu=\mu(T)$'%(T_fFD[i]))
plt.xlabel(r'$\epsilon/\epsilon_f$', fontsize=18)
plt.ylabel(r'$f_{FD}(\epsilon)$', fontsize=18)
plt.legend(fontsize=12)
plt.show()
#-----------------------------------------

# Plot of fFD with mu(T) only
#-------------------------------------------
# mu vals for the Fermi-Dirac distribution study
mu_fFD = []
for T in T_fFD:
    mu_fFD.append(fsolve(eq, x0=1, args=(T)))

for i in range(len(mu_fFD)):
    plt.plot(eps_vals, fFD(eps_vals, mu_fFD[i], T_fFD[i]), color=plot_colors[i], label=r'T=%.2f $T_F \quad \mu=\mu(T)$'%(T_fFD[i]))
plt.xlabel(r'$\epsilon/\epsilon_f$', fontsize=18)
plt.ylabel(r'$f_{FD}(\epsilon)$', fontsize=18)
plt.legend(fontsize=12)
plt.grid()
plt.show()
#-----------------------------------------