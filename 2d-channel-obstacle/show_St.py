import pandas as pd
import matplotlib.pyplot as plt

csv = pd.read_csv('data/log.strouhal')
plt.scatter(csv['Re'], csv['St'])
plt.ylim([0.13,0.2])
plt.xlabel('Re')
plt.ylabel('St')
plt.grid()
plt.show()