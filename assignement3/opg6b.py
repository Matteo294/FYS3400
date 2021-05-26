import numpy as np
import matplotlib.pyplot as plt
#from plot_set import *


def E(kx, ky):
    # notice that it should really be +- t
    E = t*np.sqrt(1 + 4*np.cos(np.sqrt(3)*ky*a/2)*np.cos(3*kx*a/2) + 4*np.cos(np.sqrt(3)*ky*a/2)**2)
    return E

def Edx(kx, ky): #dE/dk_x
    sqrt3 = np.sqrt(3)
    return t**2/(2*E(kx,ky))*(-6*a*np.cos(sqrt3*ky*a/2)*np.sin(3*kx*a/2))

def Edy(kx, ky): #dE/dk_y
    sqrt3 = np.sqrt(3)
    return t**2/(2*E(kx,ky))*(-2*sqrt3*a*np.sin(sqrt3*ky*a/2)*np.cos(3*kx*a/2) - 4*sqrt3*a*np.cos(sqrt3*ky*a/2)*np.sin(sqrt3*ky*a/2))

def m_star(kx, ky): #effective mass in [01] direction (b1,b2 basis)
    d = np.array([1/3, -1/np.sqrt(3)]) # in kx,ky basis
    return (d[0]*kx+d[1]*ky)/(d[0]*Edx(kx,ky) + d[1]*Edy(kx,ky)) #[hbar^2]

def plot3D(X,Y,Z):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    ax.set_title('surface')

def heatmap(X,Y,Z):
    plt.pcolor(X,Y,Z)
    plt.colorbar()

def plot_01_direction(): 
    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k') # [01]
    N = int(1e4+1)
    eV_to_joule = 1.602e-19
    me = 9.109e-31
    hbar = 6.582119569e-16 #[eV*s]      
    span  = 1
    k2 = np.linspace(-span, span, N)
    m = m_star(b2[0]*k2, b2[1]*k2)
    tol = 1000 # remove high peaks
    m_stripped = np.where(np.abs(m) > tol, np.nan, m) # remove high peaks
    convert = 1.60217662e-19
    m_e = 9.1093837015e-31
    m_stripped = m_stripped*convert/m_e*(1e9)**2*hbar**2
    plt.plot(k2, m_stripped, color='springgreen', linewidth=2.0)
    plt.xlabel(r"$k_2$", fontsize=16)
    plt.ylabel(r"$m^*/m_e$", fontsize=16)
    plt.grid()
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)

def m_star_gamma(): #find effective mass in limit of kx -> 0 with y = 0
    kx = 1e-20
    m = m_star(kx, 0) #[hbar^2 1/(nm^2*eV)]
    m_e = 9.1093837015e-31
    hbar = 6.582119569e-16 #[eV*s]

    #Convert
    m *= (1e9)**2   # --> [hbar^2 1/(m^2*eV)]
    m *= hbar**2    # --> [eV s**2/m**2]

    convert = 1.60217662e-19 # [eV s**2/m**2] --> [kg]
    m *= convert    # --> [kg]
    m /= m_e        # --> [m_e]
    print(m)


def plot_Bzone():
    plt.figure(num=1, dpi=80, facecolor="w", edgecolor="k") #[B-zone]
    N = int(1e2+1)

    span = 1/2*np.linalg.norm(b2)*3

    kx = np.linspace(-span,span, N) #B-zone?
    ky = np.linspace(-span,span, N) #B-zone?
    KX, KY = np.meshgrid(kx,ky)
    m3D = m_star(KX,KY)
    tol = 1000
    m3D_stipped = np.where(np.abs(m3D) > tol, np.nan, m3D) # remove high peaks

    plt.title("...")
    # plot3D(KX,KY,m3D_stipped)
    heatmap(KX, KY, m3D_stipped)
    plt.xlabel(r"$k_x$", fontsize=14)
    plt.ylabel(r"$k_y$", fontsize=14)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)

if __name__ == "__main__":
    t = -3 # [eV]
    a = 0.142 # [nm]
    # a = 0.142e-9 #[m]


    b1 = 2*np.pi/a*np.array([1/3, -1/np.sqrt(3)])
    b2 = 2*np.pi/a*np.array([1/3, -1/np.sqrt(3)])


    plot_01_direction()
    m_star_gamma()
    plot_Bzone()


    plt.show()
