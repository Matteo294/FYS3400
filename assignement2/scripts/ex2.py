import numpy as np 
from matplotlib import pyplot as plt 
from math import sin

def sinc(xarr, A=1):
    v = []
    for x in xarr:
        if x == 0:
            v.append(A*1)
        else:
            v.append(A*sin(x)/x)
    return v

xline = np.linspace(-20, 20, 10000)
plt.plot([-20, 20], [1, 1], color='gray')
plt.plot([-20, 20], [-1, -1], color='gray')
plt.plot(xline, sinc(xline, A=10) + np.cos(xline), color='black', linewidth=1.8, label=r'$\cos(x) + A \, \sin(x)/x$')
plt.tick_params(
    axis='x',      
    which='both',     
    bottom=True,   
    top=True,
    labelbottom=False,
    direction='in')
plt.tick_params(
    axis='y',      
    which='both',      
    left=True,    
    right=True,
    labelleft=False,
    direction='in')
plt.xlabel('x = ka', fontsize=14)
plt.ylabel('f(x)', fontsize=14)
plt.legend(loc=0, fontsize=14)
plt.fill_between(xline, -1, 1, facecolor='lightskyblue', alpha=0.6)
plt.show()