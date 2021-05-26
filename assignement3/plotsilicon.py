import pandas as pd 
from matplotlib import pyplot as plt

data = pd.read_csv("data.csv")
plt.semilogy(data['T'], data['n'])
plt.show()