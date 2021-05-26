import numpy as np 
from matplotlib import pyplot as plt
from numpy.core.numeric import ones
import numpy as np
import scipy.constants as scic
import scipy.integrate as scint
import scipy.optimize as scis
from scipy.integrate import *
from scipy.constants import *
from scipy.optimize import *
from sympy import symbols, Eq, solve
from sympy.functions.elementary.miscellaneous import root

eV = 1.6e-019; 
kB = 1.38e-023; # Boltzmann constant, J/mol/K
h = 6.626e-034; # Planck's constant, Js
m_0 = 9.1e-031; # Rest mass of an electron, kg
#Parameters

meff_elect = 1.08*m_0; # Effective mass of electron
meff_hole = 0.56*m_0; # Effective mass of holes 
Eg = 1.1*eV; # Energy gap of Si #setting E_v=0
Nd = 1e17
deltaE = 0.045*eV 
Ed=(Eg-deltaE) #Energy of the donor band
Nc = lambda T: 2*(2*np.pi*kB*T/h**2)**(3/2)*(meff_elect)**(3/2)/1e6

#functions
def n(Ef,T): #carrier distribution in conduction
    return ni(T)*np.exp((Ef-Efi(T))/(kB*T))
def p(Ef,T): #holes distribution in valence
    return ni(T)*np.exp((Efi(T)-Ef)/(kB*T))
def Ndp(Ef,T): # doped atoms distribution
    return Nd*(1-1/(1+0.5*np.exp((Ed-Ef)/(kB*T))))
def g(Ef,T): # neutral charge equation
    return  Ndp(Ef,T) + p(Ef,T) - n(Ef,T)
def highT(Ef,T):
    return ni(T)
def NdplowT(T):
    return np.sqrt(Nc(T)*Nd/2) * np.exp(-(Eg-Ed)/(2*kB*T))


def Efi(T):
    if isinstance(T, (int, float)):
        print("yes")
        return Eg/2 + (3/4)*kB*T*np.log(meff_hole/meff_elect)
    else:
        res = [Eg/2 + (3/4)*kB*TT*np.log(meff_hole/meff_elect) for TT in T]
        return np.array(res) #intrinsic fermi energy
ni= lambda T: 2*(2*np.pi*kB*T/h**2)**(3/2)*(meff_elect*meff_hole)**(3/4)*np.exp(-Eg/(2*kB*T)) / 10**6

Tvals = np.linspace(50,1500,1000)
Tvals1 = [1000/T for T in Tvals]
root1=np.ones(len(Tvals))
c=0
for i in Tvals:
    root1[c]=scis.fsolve(g, Eg , args=(i),xtol=eV**2)
    print(root1[c])
    c +=1 
eg= np.linspace(0,1.1,len(Tvals))
plt.plot(Tvals1,root1/eV, label=r'$E_f$', color='mediumspringgreen', linewidth=2.0)
plt.plot(Tvals1,Efi(Tvals)/eV, label=r'$E_i$', color='deepskyblue', linewidth=2.0)
plt.text(0.8, Eg/eV+0.015, "Conduction band", fontsize=10)
plt.text(0.8, Ed/eV-0.04, "Donors level", fontsize=10)
plt.text(0.8, 0.02, "Valence band", fontsize=10)
plt.plot(Tvals1,Eg*np.ones(len(Tvals))/eV,'--', color='gray')
plt.plot(Tvals1,Ed*np.ones(len(Tvals))/eV,'--', color='gray')
plt.plot(Tvals1,0*np.ones(len(Tvals)),'--', color='gray')
plt.xlabel(r"$1000/T$", fontsize=14)
plt.ylabel(r"$E \ [eV]$", fontsize=14)
plt.ylim([0, 1.3])
plt.grid()
plt.legend(fontsize=12)
plt.show()

plt.semilogy(Tvals1, n(root1,Tvals), linewidth=2.2, label='n(T)', color='black')
plt.grid()

low_T = np.linspace(50, 300, 1000)
med_T = np.linspace(150, 1500, 1000)
high_T = np.linspace(500, 1500, 1000)
lowT_x = [1000/T for T in low_T]
medT_x = [1000/T for T in med_T]
highT_x = [1000/T for T in high_T]
print(low_T)
print(Efi(low_T))

print(Ndp(Efi(low_T), low_T))
plt.semilogy(highT_x, ni(high_T), '--', label='Intrinsic', color='deepskyblue', linewidth=1.8)
plt.semilogy(medT_x, np.full(len(medT_x), Nd), '--', label='Extrinsic', color='limegreen', linewidth=1.8)
plt.semilogy(lowT_x, NdplowT(low_T), '--', label='Freeze-out', color='red', linewidth=1.8)
plt.legend(fontsize=12)
plt.ylabel(r'Concetration [$cm^{-2}$]', fontsize=14)
plt.xlabel('1000/T', fontsize=14)
plt.show()
#print(g(root1,Tvals))
#plt.semilogy(eg,n(root1,Tvals)+Ndp(root1,Tvals))



