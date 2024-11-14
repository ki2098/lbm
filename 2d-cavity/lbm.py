import numpy as np

"""
6 2 5
3 0 1
7 4 8
"""

Q = 9
D = 2

L0 = 1
N = 128
H = L0/N
C = N + 1
Dt = H
U0 = 1
Re = 1000
Nu = U0*L0/Re
Tau = 3*Nu/Dt + 0.5

E = np.array(
    [
        [ 0,  0],
        [ 1,  0],
        [ 0,  1],
        [-1,  0],
        [ 0, -1],
        [ 1,  1],
        [-1,  1],
        [-1, -1],
        [ 1, -1]
    ]
)

W = np.array([
    4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36
])

f = np.zeros((C,C,Q))
fprev = np.zeros((C,C,Q))
U = np.zeros((C,C,D))
rho = np.zeros((C,C))

rho[:,:] = 1
U[1:-1,-2,0] = U0
