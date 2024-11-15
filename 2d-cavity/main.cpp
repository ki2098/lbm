#include <cstdio>
#include <cstdlib>

const int Q = 9;
const int D = 2;
const int N = 128;
const int C = N+3;

const double L_cav_dim = 1;
const double U_cav_dim = 1;
const double Re_cav = 100;
const double nu_dim = L_cav_dim*U_cav_dim/Re_cav;
const double dx_dim = L_cav_dim/N;

const double U_cav_ndim = 0.05; // should be no more than 0.1
const double Cu_dim = U_cav_dim/U_cav_ndim;
const double dx_ndim = 1;
const double Cl_dim = dx_dim/dx_ndim;
const double Ct_dim = Cl_dim/Cu_dim;
const double dt_ndim = 1;
const double dt_dim = dt_ndim*Ct_dim;

const double L_cav_ndim = L_cav_dim/Cl_dim;
const double nu_ndim = U_cav_ndim*L_cav_ndim/Re_cav;
const double Cnu_dim = Cl_dim*Cu_dim;
const bool nu_check = (nu_ndim == nu_dim/Cnu_dim);

const double cs_ndim_sq = 1./3;
const double cs_dim_sq = cs_ndim_sq*Cu_dim*Cu_dim;

const double tau_dim  = nu_dim/cs_dim_sq + 0.5*dt_dim;
const double tau_ndim = tau_dim/Ct_dim; // should not be significantly larger than 1
const double omega = 1/tau_ndim; // better to be kept below 1.8
const bool tau_check = (nu_dim == cs_ndim_sq*(tau_ndim - 0.5*dt_dim/Ct_dim)*Cl_dim*Cl_dim/Ct_dim);

const double T_dim = 100;
const int NT = T_dim/dt_dim;

const double Re_g = U_cav_ndim*dx_ndim/nu_ndim; // should not be significantly larger than O(10)


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
double fp[C][C][Q] = {};
double U[C][C][D] = {};
double rho[C][C] = {};

double feq(double U[D], double rho, const int E[D], const double w, const double cssq) {
    double uu = U[0]*U[0] + U[1]*U[1];
    double eu = U[0]*E[0] + U[1]*E[1];
    return rho*w*(1 + eu/cssq + 0.5*(eu*eu)/(cssq*cssq) - 0.5*uu/cssq);
}

void print_info() {
    
}

void init() {
    for (int i = 1; i < C - 1; i ++) {
        U[i][C-2][0] = U_cav_ndim;
    }
    for (int i = 0; i < C; i ++) {
    for (int j = 0; j < C; j ++) {
        rho[i][j] = 1;
        for (int q = 0; q < Q; q ++) {
            f[i][j][q] = \
            feq(U[i][j], rho[i][j], E[q], W[q], cs_ndim_sq);
        }
    }}

}

void collision(
    double f[C][C][Q],
    double U[C][C][D],
    double rho[C][C],
    const int E[Q][D],
    const double W[Q],
    double cssq,
    double relax
) {
    for (int i = 0; i < C; i ++) {
    for (int j = 0; j < C; j ++) {
    for (int q = 0; q < Q; q ++) {
        f[i][j][q] += \
        (feq(U[i][j], rho[i][j], E[q], W[q], cssq) - f[i][j][q])*relax;
    }}}
}

void streaming(
    double f[C][C][Q],
    double fp[C][C][Q],
    const int E[Q][D]
) {
    for (int i = 1; i < C - 1; i ++) {
    for (int j = 1; j < C - 1; j ++) {
    for (int q = 0; q < Q; q ++) {
        f[i+E[q][0]][j+E[q][1]][q] = fp[i][j][q];
    }}}
}

void fbc(
    double f[C][C][Q],
    double u_wall
) {
    for (int j = 1; j < C - 1; j ++) {
        int I = 1;
        f[I][j][1] = f[I][j][3];
        f[I][j][5] = f[I][j][7];
        f[I][j][8] = f[I][j][6];
        I = C - 2;
        f[I][j][3] = f[I][j][1];
        f[I][j][7] = f[I][j][5];
        f[I][j][6] = f[I][j][8];
    }
    for (int i = 1; i < C - 1; i ++) {
        int J = 1;
        f[i][J][2] = f[i][J][4];
        f[i][J][6] = f[i][J][8];
        f[i][J][5] = f[i][J][7];
    }
    for (int i = 2; i < C - 2; i ++) {
        int J = C - 2;
        double r = f[i][J][0] + f[i][J][1] + f[i][J][3] + \
        2*(f[i][J][2] + f[i][J][6] + f[i][J][5]);
        f[i][J][4] = f[i][J][2];
        f[i][J][7] = f[i][J][5] + 0.5*(f[i][J][1] - f[i][J][3]) - 0.5*r*u_wall;
        f[i][J][8] = f[i][J][6] - 0.5*(f[i][J][1] - f[i][J][3]) + 0.5*r*u_wall;
    }
}

void meso_to_macro(
    double f[C][C][Q],
    double U[C][C][D],
    double rho[C][C],
    const int E[Q][D]
) {
    for (int i = 0; i < C; i ++) {
    for (int j = 0; j < C; j ++) {
        rho[i][j] = 0;
        U[i][j][0] = 0;
        U[i][j][1] = 0;
        for (int q = 0; q < Q; q ++) {
            rho[i][j] += f[i][j][q];
            U[i][j][0] += f[i][j][q]*E[q][0];
            U[i][j][1] += f[i][j][q]*E[q][1];
        }
        U[i][j][0] /= rho[i][j];
        U[i][j][1] /= rho[i][j];
    }}
}

void main_loop() {
    for (int i = 0; i < C; i ++) {
    for (int j = 0; j < C; j ++) {
    for (int q = 0; q < Q; q ++) {
        fp[i][j][q] = f[i][j][q];
    }}}

    collision(fp, U, rho, E, W, cs_ndim_sq, omega);
    streaming(f, fp, E);
    fbc(f, U_cav_ndim);
    meso_to_macro(f, U, rho, E);
}

void output() {
    FILE *file = fopen("lbm-2d-cavity.csv", "w");
    fprintf(file, "x,y,z,u,v,w,rho\n");
    for (int j = 1; j < C - 1; j ++) {
    for (int i = 1; i < C - 1; i ++) {
        fprintf(
            file,
            "%e,%e,%e,%lf,%lf,%lf,%lf\n",
            (i-1)*dx_dim, (j-1)*dx_dim, 0.0,
            U[i][j][0]*Cu_dim, U[i][j][1]*Cu_dim, 0.0,
            rho[i][j]
        );
    }}
    fclose(file);
}

int main(int argc, char **argv) {
    print_info();
    init();
    for (int t = 1; t <= NT; t ++) {
        main_loop();
        printf("\r%d/%d", t, NT);
        fflush(stdout);
    }
    printf("\n");
    output();
}