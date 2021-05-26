import numpy as np 
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

eV_to_Joule = 1.602e-19

Efact = 4.36e-18
mfact = 9.109e-31
xfact = 5.291e-11

hbar = 6.626e-34/(2*np.pi)/eV_to_Joule

m_ele = 1.08*mfact/eV_to_Joule
m_hol = 0.56*mfact/eV_to_Joule
kB = 1.381e-23/eV_to_Joule
Eg = 1.12
Ec = Eg
Ed = 45e-3
ED = Ec - Ed
Nd = 1e17
Ei0 = (Ec+ED)/2

Nc = lambda T: 2 * (m_ele*kB*T/(2*np.pi*hbar**2))**(3/2)
Nv = lambda T: 2 * (m_hol*kB*T/(2*np.pi*hbar**2))**(3/2)

def eq(x, T):
    return 2*ni(T) * np.sinh((Ei(T)-x)/(kB*T)) + Nd/(1 + 2*np.exp((x-Ei(T))/(kB*T)))

def Ei(T):
    return 0.5*Eg + 3/4*kB*T*np.log(m_hol/m_ele)

def ni(T):
    return 2* (kB*T/(2*np.pi*hbar**2))**(3/2) * (m_ele*m_hol)**(3/4) * np.exp(-Eg/(2*kB*T))

def Ef(T):
    x = fsolve(eq, 1e-3, args=(T))
    return x

def n(T):
    print(kB*T,Ef(T)-Ei(T) )
    return ni(T) * np.exp((Ef(T)-Ei(T))/(2*kB*T))

Tvals = np.linspace(100, 500, 1000)
x = []
v = []
for T in Tvals:
    x.append(1000/T)
    v.append(n(T))

plt.semilogy(x, v)
plt.show()
    