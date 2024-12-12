#include <cstdio>
#include <cstdlib>
#include <cmath>

template<typename T, int N>
struct vector_t {
    T m[N];

    T &operator[](int i) {
        return m[i];
    }
    
    const T &operator[](int i) const {
        return m[i];
    }
};

const int N = 128;
const int C = N+2;

using double9_t = vector_t<double, 9>;
using double2_t = vector_t<double, 2>;

const int VE[9][2] = {
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

const double WEIGHT[9] = {
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

const int CMK[9][2] = {
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

const double L_ = 1;
const double U_ = 1;
const double Re = 200;
const double nu_ = L_*U_/Re;
const double dx_ = L_/N;
const double dx = 1;
const double Cx_ = dx_/dx;
const double U = 0.1;
const double Cu_ = U_/U;
const double Ct_ = Cx_/Cu_;
const double dt = 1;
const double dt_ = dt*Ct_;
const double Cnu_ = Cx_*Cu_;
const double nu = nu_/Cnu_;
const double csq = 1./3.;
const double tau = nu/csq + 0.5;
const double T_ = 10;
const double NSTEP = T_/dt_;
const double Re_lattice = U*dx/nu;

void compute_statistics(
    const double9_t &f,
    double &fluxx,
    double &fluxy,
    double &density
) {
    fluxx = 0;
    fluxy = 0;
    density = 0;
    for (int fid = 0; fid < 9; fid ++) {
        fluxx += VE[fid][0]*f[fid];
        fluxy += VE[fid][1]*f[fid];
        density += f[fid];
    }
}

/**
 * c00 = rho
 * c01 = rho v
 * c10 = rho u
 * c02 = k02
 * c11 = k11
 * c20 = k20
 * c12 = k12
 * c21 = k21
 * c22 = k22 - 2(k11^2 + k20 k02)/rho
 */
void compute_cumulants(
    const double9_t *f,
    double9_t *c,
    int imax, 
    int jmax
) {
    for (int lat = 0; lat < imax*jmax; lat ++) {
        const double9_t &fc = f[lat];
        double fluxx, fluxy, density;
        compute_statistics(fc, fluxx, fluxy, density);
        double u = fluxx/density;
        double v = fluxy/density;
        double9_t kc{{}};
        for (int fid = 0; fid < 9; fid ++) {
            const int i = VE[fid][0];
            const int j = VE[fid][1];
            vector_t<double, 3> ca, cb;
            ca[0] = cb[0] = 1;
            for (int order = 1; order < 3; order ++) {
                ca[order] = ca[order - 1]*(i - u);
                cb[order] = cb[order - 1]*(j - v);
            }
            for (int kid = 0; kid < 9; kid ++) {
                const int a = CMK[kid][0];
                const int b = CMK[kid][1];
                kc[kid] += ca[a]*cb[b]*fc[fid];
            }
        }
        double9_t &cc = c[lat];
        cc[0] = density;
        cc[1] = fluxy;
        cc[3] = fluxx;
        cc[2] = kc[2];
        cc[4] = kc[4];
        cc[6] = kc[6];
        cc[5] = kc[5];
        cc[7] = kc[7];
        cc[8] = kc[8] - (2*kc[4]*kc[4] + kc[6]*kc[2])/density;
    }
}

/**
 * c00* = c00eq = c00 = rho
 * c01* = c01eq = c01 = rho v
 * c10* = c10eq = c10 = rho u
 * c02* = (1 - omega)/2 (c02 - c20) + cs^2 rho
 * c11* = (1 - omega)c11
 * c20* = (1 - omega)/2 (c20 - c02) + cs^2 rho
 * c12* = 0
 * c21* = 0
 * c22* = 0;
 */
void relax_cumulants(
    const double9_t *c,
    double9_t *cpost,
    double omega,
    int imax,
    int jmax
) {
    for (int lat = 0; lat < imax*jmax; lat ++) {
        const double9_t &cc = c[lat];
        double9_t &ccpost = cpost[lat];
        ccpost[0] = cc[0];
        ccpost[1] = cc[1];
        ccpost[3] = cc[3];
        ccpost[2] = 0.5*(1 - omega)*(cc[2] - cc[6]) + cc[0]*csq;
        ccpost[4] = (1 - omega)*cc[4];
        ccpost[6] = 0.5*(1 - omega)*(cc[6] - cc[2]) + cc[0]*csq;
        ccpost[5] = ccpost[7] = ccpost[8] = 0;
    }
}

void compute_post_pdfs(
    const double9_t *cpost,
    double9_t *fpost,
    int imax,
    int jmax
) {
    for (int lat = 0; lat < imax*jmax; lat ++) {
        const double9_t &cc = cpost[lat];
        double density = cc[0];
        double v = cc[1]/density;
        double u = cc[3]/density;
        double9_t kc{{}};
        kc[0] = density;
        kc[1] = 0;
        kc[3] = 0;
        kc[2] = cc[2];
        kc[4] = cc[4];
        kc[6] = cc[6];
        kc[5] = cc[5];
        kc[7] = cc[7];
        kc[8] = cc[8] + (2*kc[4]*kc[4] + kc[6]*kc[2])/density;

        double uu = u*u, vv = v*v, uv = u*v;
        double9_t &fc = fpost[lat];
        const double M[9][9] = {
            {  uv*(uv - u - v + 1)*0.25, u*(2*uv - u - 2*v + 1)*0.25, u*(u - 1)*0.25, v*(2*uv - 2*u - v + 1)*0.25, uv - u*0.5 - v*0.5 + 0.25, u*0.5 - 0.25, v*(v - 1)*0.25, v*0.5 - 0.25,  0.25},
            {u*(-u*vv + u + vv - 1)*0.5,                  uv*(1 - u),  u*(1 - u)*0.5,    -u*vv + u + vv*0.5 - 0.5,               v*(1 - 2*u),      0.5 - u,   0.5 - vv*0.5,           -v,  -0.5},
            {  uv*(uv + u - v - 1)*0.25, u*(2*uv + u - 2*v - 1)*0.25, u*(u - 1)*0.25, v*(2*uv + 2*u - v - 1)*0.25, uv + u*0.5 - v*0.5 - 0.25, u*0.5 - 0.25, v*(v + 1)*0.25, v*0.5 + 0.25,  0.25},
            {v*(-uu*v + uu + v - 1)*0.5,    -uu*v + uu*0.5 + v - 0.5,   0.5 - uu*0.5,                  uv*(1 - v),               u*(1 - 2*v),           -u,  v*(1 - v)*0.5,      0.5 - v,  -0.5},
            {       uu*vv - uu - vv + 1,                2*v*(uu - 1),         uu - 1,                2*u*(vv - 1),                      4*uv,          2*u,         vv - 1,          2*v,    1.},
            {v*(-uu*v - uu + v + 1)*0.5,    -uu*v - uu*0.5 + v + 0.5,   0.5 - uu*0.5,                 uv*(-v - 1),              u*(-2*v - 1),           -u, v*(-v - 1)*0.5,     -v - 0.5,  -0.5},
            {  uv*(uv - u + v - 1)*0.25, u*(2*uv - u + 2*v - 1)*0.25, u*(u + 1)*0.25, v*(2*uv - 2*u + v - 1)*0.25, uv - u*0.5 + v*0.5 - 0.25, u*0.5 + 0.25, v*(v - 1)*0.25, v*0.5 - 0.25,  0.25},
            {u*(-u*vv + u - vv + 1)*0.5,                 uv*(-u - 1), u*(-u - 1)*0.5,    -u*vv + u - vv*0.5 + 0.5,              v*(-2*u - 1),     -u - 0.5,   0.5 - vv*0.5,           -v,  -0.5},
            {  uv*(uv + u + v + 1)*0.25, u*(2*uv + u + 2*v + 1)*0.25, u*(u + 1)*0.25, v*(2*uv + 2*u + v + 1)*0.25, uv + u*0.5 + v*0.5 + 0.25, u*0.5 + 0.25, v*(v + 1)*0.25, v*0.5 + 0.25,  0.25}
        };
        for (int fid = 0; fid < 9; fid ++) {
            fc[fid] = 0;
            for (int kid = 0; kid < 9; kid ++) {
                fc[fid] += M[fid][kid]*kc[kid];
            }
        }
    }
}

void advect(
    const double9_t *fpost,
    double9_t *f,
    int imax,
    int jmax
) {
    for (int i = 1; i < imax - 1; i ++) {
    for (int j = 1; j < jmax - 1; j ++) {
    for (int q = 0; q < 9; q ++) {
        int dst = (i + VE[q][0])*jmax + (j + VE[q][1]);
        int src = i*jmax + j;
        f[dst][q] = fpost[src][q];
    }}}
}

void apply_fbc(
    double9_t *f,
    double u_wall,
    int imax,
    int jmax
) {
    /* bottom wall */
    for (int i = 1; i < imax - 1; i ++) {
        const int lat = i*jmax + 1;
        const int flist[]{2,5,8};
        for (int fid : flist) {

        }
    }
    /* left wall */
    for (int j = 1; j < jmax - 1; j ++) {
        const int lat = 1*jmax + j;
        const int flist[]{6,7,8};
        for (int fid : flist) {

        }
    }
    /* right wall */
    for (int j = 1; j < jmax - 1; j ++) {
        const int lat = (imax - 2)*jmax + j;
        const int flist[]{0,1,2};
        for (int fid : flist) {

        }
    }
}
