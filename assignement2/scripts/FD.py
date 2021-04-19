import numpy as  np 
from matplotlib import pyplot as plt 

Tf = 1
T1 = 0.01*Tf
T2 = 0.1*Tf
T3 = 0.5*Tf
T4 = 1.0*Tf
T5 = 1.2*Tf

def FD(x, T, mu=1,):
    return 1/(np.exp((x-mu)/T) +1)

xline = np.linspace(0.5, 1.5, 1000)
FD_T1 = FD(xline, T1)
FD_T2 = FD(xline, T2)
FD_T3 = FD(xline, T3)
FD_T4 = FD(xline, T4)
FD_T5 = FD(xline, T5)

plt.plot(xline, FD_T1, label=r"$0.01 \ T_f$")
plt.plot(xline, FD_T2, label=r"$0.1 \ T_f$")
plt.plot(xline, FD_T3, label=r"$0.5 \ T_f$")
plt.plot(xline, FD_T4, label=r"$1.0 \ T_f$")
plt.plot(xline, FD_T5, label=r"$1.2 \ T_f$")
plt.legend()
plt.show()