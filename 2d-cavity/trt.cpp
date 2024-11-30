#include <cstdio>
#include <cstdlib>
#include <string>

const int Q = 9;
const int D = 2;
const int N = 256;
const int C = N+2;

const double L_cav_dim = 1;
const double U_cav_dim = 1;
const double Re_cav = 5000;
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
const double tau_s_ndim = nu_ndim/cs_ndim_sq + dt_ndim*0.5; // should not be significantly larger than 1
const bool tau_s_check = (tau_s_ndim == (nu_dim/cs_dim_sq + 0.5*dt_dim)/Ct_dim);
const double Lambda = 1./4.;
const double tau_a_ndim = 0.5 + Lambda/(tau_s_ndim - 0.5);
const double T_dim = 100;
const int NT = T_dim/dt_dim;
const double Re_g = U_cav_ndim*dx_ndim/nu_ndim; // should not be significantly larger than O(10)

const double output_interval_time = 1;
const int output_interval_step = output_interval_time/dt_dim;
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

const int Link[Q] = {
    0, // 0
    3, // 1
    4, // 2
    1, // 3
    2, // 4
    7, // 5
    8, // 6
    5, // 7
    6  // 8
};

double f[C][C][Q] = {};
double fp[C][C][Q] = {};
double U[C][C][D] = {};
double rho[C][C] = {};

double feq(double U[D], double rho, const int E[D], const double w) {
    double uu = U[0]*U[0] + U[1]*U[1];
    double eu = U[0]*E[0] + U[1]*E[1];
    return rho*w*(1 + eu/cs_ndim_sq + 0.5*(eu*eu)/(cs_ndim_sq*cs_ndim_sq) - 0.5*uu/cs_ndim_sq);
}

void print_info() {
    
}

void init() {
    // for (int i = 1; i < C - 1; i ++) {
    //     U[i][C-2][0] = U_cav_ndim;
    // }
    for (int i = 0; i < C; i ++) {
    for (int j = 0; j < C; j ++) {
        U[i][j][0] = 0;
        U[i][j][1] = 0;
        rho[i][j] = 1;
        for (int q = 0; q < Q; q ++) {
            f[i][j][q] = feq(U[i][j], rho[i][j], E[q], W[q]);
            fp[i][j][q] = f[i][j][q];
        }
    }}
    #pragma acc enter data \
    copyin(E, W, f, fp, U, rho, Link)
}

void finalize() {
    #pragma acc exit data \
    delete(E, W, f, fp, U, rho, Link)
}

// TRT collision operator
void collision(
    double f[C][C][Q],
    double fp[C][C][Q],
    double U[C][C][D],
    double rho[C][C],
    double tau_s,
    double tau_a
) {
    #pragma acc parallel loop independent collapse(3) \
    present(f, fp, U, rho, E, W, Link) \
    firstprivate(tau_s, tau_a)
    for (int i = 0; i < C; i ++) {
    for (int j = 0; j < C; j ++) {
    for (int q = 0; q < Q; q ++) {
        int l = Link[q];
        double f_s = 0.5*(f[i][j][q] + f[i][j][l]);
        double f_a = 0.5*(f[i][j][q] - f[i][j][l]);
        double feqq = feq(U[i][j], rho[i][j], E[q], W[q]);
        double feql = feq(U[i][j], rho[i][j], E[l], W[l]);
        double feq_s = 0.5*(feqq + feql);
        double feq_a = 0.5*(feqq - feql);
        fp[i][j][q] = f[i][j][q] + (feq_s - f_s)/tau_s + (feq_a - f_a)/tau_a;
    }}}
}

void streaming(
    double f[C][C][Q],
    double fp[C][C][Q]
) {
    #pragma acc parallel loop independent collapse(3) \
    present(f, fp, E)
    for (int i = 1; i < C - 1; i ++) {
    for (int j = 1; j < C - 1; j ++) {
    for (int q = 0; q < Q; q ++) {
        f[i+E[q][0]][j+E[q][1]][q] = fp[i][j][q];
    }}}
}

void fbc(
    double f[C][C][Q],
    double rho[C][C],
    double u_wall
) {
    #pragma acc parallel loop independent \
    present(f, rho, E, W) \
    firstprivate(u_wall)
    for (int i = 1; i < C - 1; i ++) {
        const int j = C - 2; // top lid
        f[i][j][7] = f[i+1][j+1][5] + E[7][0]*2*rho[i][j]*u_wall*W[7]/cs_ndim_sq;
        f[i][j][4] = f[i  ][j+1][2];
        f[i][j][8] = f[i-1][j+1][6] + E[8][0]*2*rho[i][j]*u_wall*W[8]/cs_ndim_sq;
    }
    #pragma acc parallel loop independent \
    present(f)
    for (int i = 1; i < C - 1; i ++) {
        const int j = 1; // bottom wall
        f[i][j][6] = f[i+1][j-1][8];
        f[i][j][2] = f[i  ][j-1][4];
        f[i][j][5] = f[i-1][j-1][7];
    }
    #pragma acc parallel loop independent \
    present(f)
    for (int j = 1; j < C - 1; j ++) {
        const int i = C - 2; // right wall
        f[i][j][6] = f[i+1][j-1][8];
        f[i][j][3] = f[i+1][j  ][1];
        f[i][j][7] = f[i+1][j+1][5];
    }
    #pragma acc parallel loop independent \
    present(f)
    for (int j = 1; j < C - 1; j ++) {
        const int i = 1; // left wall
        f[i][j][5] = f[i-1][j-1][7];
        f[i][j][1] = f[i-1][j  ][3];
        f[i][j][8] = f[i-1][j+1][6];
    }
}

void meso_to_macro(
    double f[C][C][Q],
    double U[C][C][D],
    double rho[C][C],
    const int E[Q][D]
) {
    #pragma acc parallel loop independent collapse(2) \
    present(f, U, rho, E)
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
    // #pragma acc parallel loop independent collapse(3) \
    // present(fp, f)
    // for (int i = 0; i < C; i ++) {
    // for (int j = 0; j < C; j ++) {
    // for (int q = 0; q < Q; q ++) {
    //     fp[i][j][q] = f[i][j][q];
    // }}}

    collision(f, fp, U, rho, tau_s_ndim, tau_a_ndim);
    streaming(f, fp);
    fbc(f, rho, U_cav_ndim);
    meso_to_macro(f, U, rho, E);
}

void output(int n) {
    #pragma acc update \
    self(U, rho)
    std::string fname = "data/trt.csv." + std::to_string(n);
    FILE *file = fopen(fname.c_str(), "w");
    fprintf(file, "x,y,z,u,v,w,rho\n");
    for (int j = 1; j < C - 1; j ++) {
    for (int i = 1; i < C - 1; i ++) {
        fprintf(
            file,
            "%e,%e,%e,%lf,%lf,%lf,%lf\n",
            (i-1+0.5)*dx_dim, (j-1+0.5)*dx_dim, 0.0,
            U[i][j][0]*Cu_dim, U[i][j][1]*Cu_dim, 0.0,
            rho[i][j]
        );
    }}
    fclose(file);
}

void output() {
    #pragma acc update \
    self(U, rho)
    FILE *file = fopen("data/trt.csv", "w");
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
    for (int step = 1; step <= NT; step ++) {
        main_loop();
        printf("\r%d/%d", step, NT);
        fflush(stdout);
        if (step % output_interval_step == 0) {
            output(step / output_interval_step);
        }
    }
    // output();
    printf("\n");
}