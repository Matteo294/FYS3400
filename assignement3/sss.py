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
h= 6.626e-034; # Planck's constant, Js
m_0 = 9.1e-031; # Rest mass of an electron, kg
#Parameters

meff_elect = 1.08*m_0; # Effective mass of electron
meff_hole = 0.56*m_0; # Effective mass of holes 
Eg = 1.1*eV; # Energy gap of Si #setting E_v=0
Nd = 1e17
deltaE = 0.045*eV 
Ed=(Eg-deltaE) #Energy of the donor band


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


Efi= lambda T: Eg/2 + (3/4)*kB*T*np.log(meff_hole/meff_elect) #intrinsic fermi energy
ni= lambda T: 2*(2*np.pi*kB*T/h**2)**(3/2)*(meff_elect*meff_hole)**(3/4)*np.exp(-Eg/(2*kB*T)) / 10**6

Tvals = np.linspace(70,1300,1000)
Tvals1 = [1000/T for T in Tvals]
root1=np.ones(len(Tvals))
c=0
for i in Tvals:
    root1[c]=scis.fsolve(g, Eg , args=(i),xtol=eV**2)
    print(root1[c])
    c +=1 
eg= np.linspace(0,1.1,len(Tvals))
plt.plot(Tvals1,root1/eV)
plt.plot(Tvals1,Efi(Tvals)/eV)
plt.plot(Tvals1,Eg*np.ones(len(Tvals))/eV,'--')
plt.plot(Tvals1,Ed*np.ones(len(Tvals))/eV,'--')
plt.plot(Tvals1,0*np.ones(len(Tvals)),'--')
plt.xlabel(r"$1000/T$")
plt.ylabel(r"$E [eV]$")
plt.grid()
plt.show()
plt.semilogy(Tvals1, n(root1,Tvals))
plt.grid()

low_T = np.linspace(100, 300, 1000)
med_T = np.linspace(300, 500, 1000)
high_T = np.linspace(500, 1500, 1000)
low_x = [1000/T for T in low_T]
med_x = [1000/T for T in med_T]
high_x = [1000/T for T in high_T]

plt.plot(low_x, Ndp(E))

plt.show()
#print(g(root1,Tvals))
#plt.semilogy(eg,n(root1,Tvals)+Ndp(root1,Tvals))



