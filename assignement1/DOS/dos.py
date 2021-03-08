import numpy as np 
from matplotlib import pyplot as plt 
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

w = np.linspace(0, 0.95, 10000)
plt.plot(w, np.full(10000, 1.5), linewidth=1.8, color='red', label='DOS(k)')
plt.plot(w, 1/np.sqrt(1-w**2), linewidth=1.8, color='deepskyblue', label=r'DOS($\omega$)')
ax = plt.gca()
ax.axes.xaxis.set_ticklabels([], fontsize=14)
ax.axes.yaxis.set_ticklabels([], fontsize=14)
plt.title('Density of states', fontsize=14)
plt.legend()
plt.show()
