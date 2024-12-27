import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

csv = pd.read_csv('data/cd_history.csv')
avg = np.average(csv['cd'])
std = np.std(csv['cd'])
print('avg=%s std=%s'%(avg, std))

plt.plot(csv['t'], csv['cd'])
plt.show()
