import pandas as pd
import matplotlib.pyplot as plt

csv = pd.read_csv('data/log.StCD')

plt.figure()
plt.scatter(csv['Re'], csv['St'], marker='x', c='red')
plt.xlim([0, 450])
plt.ylim([0.1, 0.17])
plt.grid()
plt.xlabel('Re')
plt.ylabel('St')
plt.show()