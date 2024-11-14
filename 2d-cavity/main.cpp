#include <cstdio>
#include <cstdlib>

const int Q = 9;
const int D = 2;
const double L0 = 1;
const double U0 = 1;
const double Re = 1000;
const double Nu = U0*L0/Re;
const int N  = 128;
const int C = N + 1;
const double Delta = L0/N; 
const double Dt = Delta;
const double Tau = 3*Nu/Dt + 0.5;

/**
 * 6 2 5
 * 3 0 1
 * 7 4 8
 */

const int E[Q][D] = {
    { 0,  0},
    { 1,  0},
    { 0,  1},
    {-1,  0},
    { 0, -1},
    { 1,  1},
    {-1,  1},
    {-1, -1},
    { 1, -1}
};

const double W[Q] = {
    4./9,
    1./9,
    1./9,
    1./9,
    1./9,
    1./36,
    1./36,
    1./36,
    1./36
};

double f[C][C][Q] = {};
double ftmp[C][C][Q] = {};
double U[C][C][D] = {};
double rho[C][C] = {};

double feq(double U[D], double rho, double E[D], double w) {
    double uu = U[0]*U[0] + U[1]*U[1];
    double eu = U[0]*E[0] + U[1]*E[1];
    return rho*w*(1 + 3*eu + 4.5*eu*eu - 1.5*uu);
}

void collision(
    double f[C][C][Q],
    double U[C][C][D],
    double rho[C][C],
    double E[Q][D],
    double W[Q],
    double tau
) {
    for (int i = 0; i < C; i ++) {
    for (int j = 0; j < C; j ++) {
    for (int q = 0; q < Q; q ++) {
        f[i][j][q] += \
        (feq(U[i][j], rho[i][j], E[q], W[q]) - f[i][j][q])/tau;
    }}}
}