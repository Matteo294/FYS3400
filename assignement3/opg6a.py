import numpy as np
import matplotlib.pyplot as plt 
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def reciprocal_2D(a):
    b =(2*np.pi/a_0)*a.I
    return b 
def E(kx,ky):
    # notice that it should really be +- t
    E = t*np.sqrt(1 + 4*np.cos(np.sqrt(3)*ky*a_0/2)*np.cos(3*kx*a_0/2) + 4*np.cos(np.sqrt(3)*ky*a_0/2)**2)
    return E

def Edx(kx, ky): #dE/dk_x
    sqrt3 = np.sqrt(3)
    return t**2/(2*E(kx,ky))*(-6*a*np.cos(sqrt3*ky*a/2)*np.sin(3*kx*a/2))

def Edy(kx, ky): #dE/dk_y
    sqrt3 = np.sqrt(3)
    return t**2/(2*E(kx,ky))*(-2*sqrt3*a*np.sin(sqrt3*ky*a/2)*np.cos(3*kx*a/2) - 4*sqrt3*a*np.cos(sqrt3*ky*a/2)*np.sin(sqrt3*ky*a/2))

def mass(x, kx, ky, vx, vy):
    return x * hbar**2 * (vx*kx + vy*ky) / (vx*Edx + vy*Edy) 

###
#define the basis
a_0=0.142e-9
hbar=6.626e-34
t=3
a1=[0.5*3*a_0,0.5*np.sqrt(3)*a_0]
a2=[0.5*3*a_0,-0.5*np.sqrt(3)*a_0]
A=np.matrix([a1,a2])
B=reciprocal_2D(A)
#b1=np.squeeze(np.asarray(B[:,0]))
#b2=np.squeeze(np.asarray(B[:,1]))
b1=[(2*np.pi/a_0)/3,-(2*np.pi/a_0)/np.sqrt(3)]
b2=[(2*np.pi/a_0)/3,(2*np.pi/a_0)/np.sqrt(3)]

### plot
fig, axs = plt.subplots(3, sharex=True, sharey=True)
span=1#np.pi/a_0
x=np.linspace(-span,span,10000)
axs[0].set_title(r'Direction $\mathbf{b}_2$')
y1=E(b2[0]*x,b2[1]*x)
axs[0].plot(x,y1, color='mediumspringgreen', linewidth=2.0, label='Conduction band')
axs[0].plot(x,-y1, color='royalblue', linewidth=2.0, label='Valence band')
y2=E(b1[0]*x,b1[1]*x)
axs[1].set_title(r'Direction $\mathbf{b}_1$')
axs[1].plot(x,y2, color='mediumspringgreen', linewidth=2.0, label='Conduction band')
axs[1].plot(x,-y2, color='royalblue', linewidth=2.0, label='Valence band')
y3=E((b1[0]+b2[0])*x,(b1[1]+b2[1])*x)
print(y3)
axs[2].set_title(r'Direction $\mathbf{b}_1 + \mathbf{b}_2$')
axs[2].plot(x,y3, color='mediumspringgreen', linewidth=2.0, label='Conduction band')
axs[2].plot(x,-y3, color='royalblue', linewidth=2.0, label='Valence band')
plt.setp(axs[:], ylabel=r'E($\mathbf{k}$)')
plt.setp(axs[2], xlabel='k')
#grid
axs[0].grid()
axs[1].grid()
axs[2].grid()
for i in range (0,3):
    axs[i].plot([0,0],[-E(0,0), +E(0,0)], '--', label='Max gap = 18 eV', color='red')
    axs[i].plot([0.5,0.5],[-E(b2[0]*0.5,b2[1]*0.5), +E(b2[0]*0.5,b2[1]*0.5)], '--', label='Min gap = 6 eV', color='orange')
    axs[i].legend(loc='center right')
#axs[1].plot([-0.5 -0.5],)
#axs[2].plot([-0.5 -0.5],)

#axs[2].xlabel(r"$k_2$", fontsize=14)
#axs[0].ylabel(r"$m^*$", fontsize=14)
plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
plt.show()

xline = np.linspace(-10, 10, 1000)
plt.plot(xline, mass(xline, ))