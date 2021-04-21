import numpy as  np 
from matplotlib import pyplot as plt 
from scipy.optimize import fsolve

def DOS(E, beta):
    return np.log(E**3/np.tanh(beta*E/2))

def Cv(E, beta):
    return (beta*E)**3 * np.exp(beta*E)/(np.exp(beta*E) + 1)**2

temperatures = [1e-2]

for T in temperatures:
    beta = 1/T 
    Eline = np.linspace(1e-10, 15*T, 10000)
    val = Cv(Eline, beta)

    maxVal = max(val)
    idxMax = np.where(val == max(val))[0][0]
    maxVal_point = Eline[idxMax]

    diff = val - maxVal/5

    idx_FWHM_point = np.where(abs(diff) == min(abs(diff)))[0][0]
    FWHM_point = Eline[idx_FWHM_point]

    
    diff2 = np.delete(diff, range(idx_FWHM_point-5,idx_FWHM_point+5))
    idx_FWHM_point_2 = np.where(abs(diff2) == min(abs(diff2)))[0][0]
    FWHM_point_2 = Eline[idx_FWHM_point_2]

    plt.plot(Eline, val, color='black', linewidth=2.0, label=r'$c_v(\epsilon$)')
    plt.plot([maxVal_point, maxVal_point], [0, maxVal], '--', linewidth=0.8, color='black', alpha=0.5)
    plt.plot([FWHM_point, FWHM_point], [0, Cv(FWHM_point, beta)], '--', linewidth=0.8, color='black', alpha=0.5)
    plt.plot([FWHM_point_2, FWHM_point_2], [0, Cv(FWHM_point_2, beta)], '--', linewidth=0.8, color='black', alpha=0.5)
    plt.fill_between(Eline[idx_FWHM_point_2:idx_FWHM_point], 0, val[idx_FWHM_point_2:idx_FWHM_point], facecolor='black', alpha=0.3)
    plt.plot([1e-10, maxVal_point], [maxVal, maxVal], '--', color='black', linewidth=0.8)
    plt.plot([1e-10, FWHM_point_2], [Cv(FWHM_point_2, beta), Cv(FWHM_point_2, beta)], '--', color='black', linewidth=0.8)
    plt.legend(fontsize=14)
    plt.xlabel(r'$\epsilon$', fontsize=18)
    plt.ylabel(r'$c_v(\epsilon$)', fontsize=18)
    plt.show()
    
    plt.plot(Eline, DOS(Eline, beta), color='black', linewidth=2.0, label=r'DOS($\epsilon$)')
    plt.plot([maxVal_point, maxVal_point], [min(DOS(Eline, beta)), DOS(maxVal_point, beta)], '--', color='black', linewidth=0.8)
    plt.plot([FWHM_point, FWHM_point], [min(DOS(Eline, beta)), DOS(FWHM_point, beta)], '--', color='black', linewidth=0.8)
    plt.plot([FWHM_point_2, FWHM_point_2], [min(DOS(Eline, beta)), DOS(FWHM_point_2, beta)], '--', color='black', linewidth=0.8)
    plt.fill_between(Eline[idx_FWHM_point_2:idx_FWHM_point], min(DOS(Eline, beta)), DOS(Eline[idx_FWHM_point_2:idx_FWHM_point], beta), facecolor='black', alpha=0.3)
    plt.plot([1e-10, FWHM_point_2], [DOS(FWHM_point_2, beta), DOS(FWHM_point_2, beta)], '--', color='black', linewidth=0.8)
    plt.plot([1e-10, FWHM_point], [DOS(FWHM_point, beta), DOS(FWHM_point, beta)], '--', color='black', linewidth=0.8)
    plt.plot([1e-10, maxVal_point], [DOS(maxVal_point, beta), DOS(maxVal_point, beta)], '--', color='black', linewidth=0.8)
    plt.legend(fontsize=14)
    plt.xlabel(r'$\epsilon$', fontsize=18)
    plt.ylabel(r'DOS($\epsilon$)', fontsize=18)
    plt.show()