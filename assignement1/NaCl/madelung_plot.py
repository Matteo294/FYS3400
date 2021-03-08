from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

data = pd.read_csv('madelung.csv')
plt.plot(data['N_atoms'], data['Madelung'], linewidth=2.2, color='cornflowerblue', label='Model 1')
plt.plot(data['N_atoms'], data['Madelung_borders'], linewidth=2.2, color='limegreen', label='Model 2')
plt.plot(data['N_atoms'], np.full(len(data['N_atoms']), 2*np.log(2)), linewidth=2.2, color='black', label=r'$2 \, \log 2$')
plt.xlabel('Number of atoms', fontsize=16)
plt.ylabel('Madelung constant', fontsize=16)
plt.legend()
plt.grid()
plt.show()