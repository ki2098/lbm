import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv = pd.read_csv('data/probe.csv')

value = np.asarray(csv['value'])
minimum = min(value)
maximum = max(value)
center = 0.5*(minimum + maximum)
value = value - center
amp = maximum - minimum
value = value / amp

print(minimum, maximum)

fft = np.fft.fft(value)
n = len(fft)
fft = abs(fft[0:n//2])

fft_sorted = np.sort(fft)
print('k\ts')
for i in range(-1,-11,-1):
    k = np.where(fft == fft_sorted[i])[0][0]
    s = k/100
    print('%s\t%s'%(k,s))


