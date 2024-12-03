#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <Eigen/Eigen>

using namespace std;

const int Q = 9;
const int D = 2;
const int N = 128;
const int C = N+2;

const double L_cav_dim = 1;
const double U_cav_dim = 1;
const double Re_cav = 3200;
const double nu_dim = L_cav_dim*U_cav_dim/Re_cav;
const double dx_dim = L_cav_dim/N;
const double U_cav_ndim = 0.1; // should be no more than 0.1
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

const double output_interval_time = 1;
const int output_interval_step = output_interval_time/dt_dim;

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

const double s0 = 1, s3 = 1, s5 = 1;
const double s7 = 1/tau_ndim, s8 = 1/tau_ndim;
const double s4 = 1.2, s6 = 1.2;
const double s1 = 1.1;
const double s2 = 1;

double f[C][C][Q] = {};
double f_new[C][C][Q] = {};
double m[C][C][Q] = {};
double m_new[C][C][Q] = {};
double U[C][C][D] = {};
double rho[C][C] = {};
double M[Q][Q];
double M_inv[Q][Q];
double S_hat[Q] = {s0, s1, s2, s3, s4, s5, s6, s7, s8};

void init_matrices() {
    Eigen::Matrix<double, Q, Q> _M {
        { 1,  1,  1,  1,  1,  1,  1,  1,  1},
        {-4, -1, -1, -1, -1,  2,  2,  2,  2},
        { 4, -2, -2, -2, -2,  1,  1,  1,  1},
        { 0,  1,  0, -1,  0,  1, -1, -1,  1},
        { 0, -2,  0,  2,  0,  1, -1, -1,  1},
        { 0,  0,  1,  0, -1,  1,  1, -1, -1},
        { 0,  0, -2,  0,  2,  1,  1, -1, -1},
        { 0,  1, -1,  1, -1,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  1, -1,  1, -1}
    };
    Eigen::Matrix<double, Q, Q> _M_inv = _M.inverse();
    Eigen::DiagonalMatrix<double, Q, Q> _S_hat{s0, s1, s2, s3, s4, s5, s6, s7, s8};
    cout << "M = \n" << _M << endl;
    cout << "M inverse = \n" << _M_inv << endl;
    cout << "S hat = \n" << static_cast<Eigen::Matrix<double, Q, Q>>(_S_hat) << endl;
    Eigen::Matrix<double, Q, Q> _S = _M_inv*_S_hat*_M;
    cout << "S = \n" << _S << endl;
    Eigen::EigenSolver<Eigen::Matrix<double, Q, Q>> _ei(_S);
    cout << "S eigen = \n" << _ei.eigenvalues() << endl;

    for (int i = 0; i < Q; i ++) {
    for (int j = 0; j < Q; j ++) {
        M[i][j] = _M(i, j);
        M_inv[i][j] = _M_inv(i, j);
    }}
}

void apply_matrix(double dst[C][C][Q], double Mat[Q][Q], double src[C][C][Q]) {
    for (int i = 0; i < C; i ++) {
    for (int j = 0; j < C; j ++) {
    for (int q = 0; q < Q; q ++) {
        double sum = 0;
        for (int k = 0; k < Q; k ++) {
            sum += Mat[q][k]*src[i][j][k];
        }
        dst[i][j][q] = sum;
    }}}
}

void get_meq(double meq[Q], double U[D], double rho) {
    double u = U[0];
    double v = U[1];
    meq[0] = rho;
    meq[1] = rho*(- 2 + 3*(u*u + v*v));
    meq[2] = rho*(  1 - 3*(u*u + v*v));
    meq[3] = rho*u;
    meq[4] = rho*(- u);
    meq[5] = rho*v;
    meq[6] = rho*(- v);
    meq[7] = rho*(u*u - v*v);
    meq[8] = rho*u*v;
}

void collide(double m_new[C][C][Q], double m[C][C][Q], double S_hat[Q], double U[C][C][D], double rho[C][C]) {
    for (int i = 0; i < C; i ++) {
    for (int j = 0; j < C; j ++) {
        double meq[9];
        get_meq(meq, U[i][j], rho[i][j]);
        for (int q = 0; q < Q; q ++) {
            m_new[i][j][q] = m[i][j][q] + S_hat[q]*(meq[q] - m[i][j][q]);
        }
    }}
}

void stream(double f_new[C][C][Q], double f[C][C][Q]) {
    for (int i = 1; i < C - 1; i ++) {
    for (int j = 1; j < C - 1; j ++) {
    for (int q = 0; q < Q; q ++) {
        f_new[i + E[q][0]][j + E[q][1]][q] = f[i][j][q];
    }}}
}

void apply_fbc(double f[C][C][Q], double rho[C][C], double u_wall) {
    for (int i = 1; i < C - 1; i ++) {
        const int j = 1; // bottom wall
        f[i][j][6] = f[i+1][j-1][8];
        f[i][j][2] = f[i  ][j-1][4];
        f[i][j][5] = f[i-1][j-1][7];
    }

    for (int j = 1; j < C - 1; j ++) {
        const int i = C - 2; // right wall
        f[i][j][6] = f[i+1][j-1][8];
        f[i][j][3] = f[i+1][j  ][1];
        f[i][j][7] = f[i+1][j+1][5];
    }

    for (int j = 1; j < C - 1; j ++) {
        const int i = 1; // left wall
        f[i][j][5] = f[i-1][j-1][7];
        f[i][j][1] = f[i-1][j  ][3];
        f[i][j][8] = f[i-1][j+1][6];
    }

    for (int i = 1; i < C - 1; i ++) {
        const int j = C - 2; // top lid
        f[i][j][7] = f[i+1][j+1][5] + E[7][0]*2*rho[i][j]*u_wall*W[7]/cs_ndim_sq;
        f[i][j][4] = f[i  ][j+1][2];
        f[i][j][8] = f[i-1][j+1][6] + E[8][0]*2*rho[i][j]*u_wall*W[8]/cs_ndim_sq;
    }
}

void get_macro(double f[C][C][Q], double U[C][C][D], double rho[C][C]) {
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
    apply_matrix(m, M, f);
    collide(m_new, m, S_hat, U, rho);
    apply_matrix(f_new, M_inv, m_new);
    stream(f, f_new);
    apply_fbc(f, rho, U_cav_ndim);
    get_macro(f, U, rho);
}

void init() {
    for (int i = 0; i < C; i ++) {
    for (int j = 0; j < C; j ++) {
        U[i][j][0] = 0;
        U[i][j][1] = 0;
        rho[i][j] = 1;
        get_meq(m[i][j], U[i][j], rho[i][j]);
    }}
    init_matrices();
    apply_matrix(f, M_inv, m);
    get_macro(f, U, rho);
}

void output() {
    FILE *file = fopen("data/mrt.csv", "w");
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
    init();
    for (int step = 1; step <= NT; step ++) {
        main_loop();
        printf("\r%d/%d", step, NT);
        fflush(stdout);
        if (step % output_interval_step == 0) {
            // output(step / output_interval_step);
        }
    }
    output();
    printf("\n");
    return 0;
}
