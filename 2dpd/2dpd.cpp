#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <chrono>
#include <iostream>
#include <fstream>

#include "../vector_type.h"

using double9_t = vector_t<double, 9>;
using double2_t = vector_t<double, 2>;
using double3_t = vector_t<double, 3>;

const int Ve[9][2] = {
    {-1, -1},
    {-1,  0},
    {-1,  1},
    { 0, -1},
    { 0,  0},
    { 0,  1},
    { 1, -1},
    { 1,  0},
    { 1,  1}
};

const double Wgt[9] = {
    1./36,
    1./9,
    1./36,
    1./9,
    4./9,
    1./9,
    1./36,
    1./9,
    1./36
};

const int Cmk[9][2] = {
    {0, 0},
    {0, 1},
    {0, 2},
    {1, 0},
    {1, 1},
    {1, 2},
    {2, 0},
    {2, 1},
    {2, 2}
};
#pragma acc declare \
copyin(Ve, Wgt, Cmk)

const double L_ = 1;
const double U_ = 1;
const double Re = 50;
const double Lx_ = 50*L_;
const double Ly_ = 25*L_;
const int Cells_per_length = 20;
const int Ghost_cell = 1;
const double nu_ = L_*U_/Re;
const double dx_ = 1./Cells_per_length;
const double dx = 1;
const double Cx_ = dx_/dx;
const double U = 0.05;
const double Cu_ = U_/U;
const double Ct_ = Cx_/Cu_;
const double dt = 1;
const double dt_ = dt*Ct_;
const double Cnu_ = Cx_*Cu_;
const double nu = nu_/Cnu_;
const double csq = 1./3.;
const double csqi = 3;
const double tau = nu*csqi + 0.5;
const double omega = 1/tau;
const double Re_lattice = U*dx/nu;
const double CFL = U_*dt_/dx_;
const double Cpressure_ = 1.*Cu_*Cu_;
const double CDD = 13;

class Cumu {
public:
    double9_t *f, *fpost, *c, *cpost;
    double2_t *shift;
    double3_t *mac; /* [ρ u v] */

    int imax, jmax;
    double omega;

    Cumu(int imax, int jmax, double omega) : imax(imax), jmax(jmax), omega(omega) {
        f = new double9_t[imax*jmax];
        fpost = new double9_t[imax*jmax];
        c = new double9_t[imax*jmax];
        cpost = new double9_t[imax*jmax];
        shift = new double2_t[imax*jmax];
        mac = new double3_t[imax*jmax];

        #pragma acc enter data \
        copyin(this[0:1], f[:imax*jmax], fpost[:imax*jmax], c[:imax*jmax], cpost[:imax*jmax], shift[:imax*jmax], mac[:imax*jmax])
    }

    ~Cumu() {
        delete[] f;
        delete[] fpost;
        delete[] c;
        delete[] cpost;
        delete[] shift;
        delete[] mac;

        #pragma acc exit data \
        delete(f[:imax*jmax], fpost[:imax*jmax], c[:imax*jmax], cpost[:imax*jmax], shift[:imax*jmax], mac[:imax*jmax], this[0:1])
    }

    void print_info() {
        printf("CUMULANT LBM\n");
        printf("\tdomain size = (%d %d)\n", imax, jmax);
        printf("\tghost cell = %d\n", Ghost_cell);
        printf("\trelaxation rate = %lf\n", omega);
        printf("\tRe = %lf\n", Re);
    }
};

double3_t get_macroscopic_val(const double9_t &f, const double2_t &shift) {
    double density = 0, u = 0, v = 0;
    for (int fid = 0; fid < 9; fid ++) {
        density += f[fid];
        u += Ve[fid][0]*f[fid];
        v += Ve[fid][1]*f[fid];
    }
    u = (u + shift[0])/density;
    v = (v + shift[1])/density;
    return double3_t{{density, u, v}};
}

double9_t get_equilibrium_cumulant(const double3_t &mac) {
    double density = mac[0];
    double u = mac[1];
    double v = mac[2];
    return double9_t{{
        /* 00       01          02          10     11 12      20      21 22 */
        density, density*v, csq*density, density*u, 0, 0, csq*density, 0, 0
    }};
}
