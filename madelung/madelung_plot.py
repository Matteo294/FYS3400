from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

data = pd.read_csv('madelung.csv')
plt.plot(data['N_molecules'], data['Madelung'], linewidth=1.8, color='royalblue', label='N atoms')
plt.plot(data['N_molecules'], np.full(len(data['N_molecules']), 2*np.log(2)), linewidth=2.2, color='black', label=r'$2 \, \log 2$')
plt.plot(data['N_molecules'], data['Madelung_borders'], linewidth=1.8, color='green', label='N atoms')
plt.xlabel('Number of molecules', fontsize=16)
plt.ylabel('Madelung constant', fontsize=16)
plt.legend()
plt.title('Madelung constant as a function of the number of molecules in the lattice', fontsize=18)
plt.show()