import numpy as np
import numpy.typing as npt
import math

Q = 9
D = 2
L0 = 1
U0 = 1
Re = 1000
Nu = U0*L0/Re
N  = 128
Delta = L0/N 
Dt = Delta
Tau = 3*Nu/Dt + 0.5
C = N + 1

def info():
    print('Sim info:')
    print(f'\tmodel = D{D}Q{Q}')
    print(f'\tL0 = {L0}')
    print(f'\tU0 = {U0}')
    print(f'\tRe = {Re}')
    print(f'\tNu = {Nu}')
    print(f'\tnon-dimensional Tau = {Tau}')
    print(f'\tDx = Dy = {Delta}')
    print(f'\tDt = {Dt}')
    print(f'\tgrid size = ({C} {C})')

# 6 2 5
# 3 0 1
# 7 4 8
W = np.array([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36])
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

f = np.zeros((C, C, Q))
ftmp = np.zeros((C, C, Q))
U = np.zeros((C, C, D))
rho = np.zeros((C, C))

def feq(U, rho, e, w):
    uu = U[0]**2 + U[1]**2
    eu = U[0]*e[0] + U[1]*e[1]
    return rho*w*(1 + 3*eu + 4.5*eu**2 - 1.5*uu)

def collision(f, U, rho, E, W):
    for i in range(C):
        for j in range(C):
            for q in range(Q):
                f[i, j, q] += (feq(W[q], E[q, :], U[i, j, :], rho[i, j]) - f[i, j, q])/Tau


if __name__ == "__main__":
    info()
