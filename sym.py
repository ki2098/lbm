from sympy import *

a, b, u, v, i, j = symbols('a b u v i j')

expr = ((i-u)**a)*((j-v)**b)

mat = zeros(9)

AB = [
    [0,0],
    [1,0],
    [0,1],
    [2,0],
    [1,1],
    [0,2],
    [2,1],
    [1,2],
    [2,2]
]

IJ = [
    [ 0, 0],
    [ 1, 0],
    [ 0, 1],
    [-1, 0],
    [ 0,-1],
    [ 1, 1],
    [-1, 1],
    [-1,-1],
    [ 1,-1]
]

for kid in range(9):
    for fid in range(9):
        mat[kid,fid] = expr.subs([(a, AB[kid][0]), (b, AB[kid][1]), (i, IJ[fid][0]), (j, IJ[fid][1])])

mat_i = mat.inv()

mat_t = expand(mat*mat_i)
