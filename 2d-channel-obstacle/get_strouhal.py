import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv = pd.read_csv('data/probe.csv')

value = np.asarray(csv['value'])
minimum = min(value)
maximum = max(value)
amp = (maximum - minimum)/2
value = (value - (maximum + minimum)/2)/amp
fft = 2*np.fft.fft(value)/value.size
power = (np.abs(fft))**2
dt = 0.0025
freq = np.fft.fftfreq(value.size, dt)
# idx = np.argsort(freq)

fig, ax = plt.subplots(2, 1)
ax[0].plot(csv['t'], value)
ax[1].plot(freq[:500], power[:500])
ax[1].grid()
plt.show()


