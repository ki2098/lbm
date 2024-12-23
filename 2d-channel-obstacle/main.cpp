#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <chrono>

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

using double9_t = vector_t<double, 9>;
using double2_t = vector_t<double, 2>;

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
const double Re = 5.8e5;
const double Lx_ = 20*L_;
const double Ly_ = 5*L_;
const int Cells_per_length = 10;
const int Ghost_cell = 1;
const double nu_ = L_*U_/Re;
const double dx_ = 1./Cells_per_length;
const double dx = 1;
const double Cx_ = dx_/dx;
const double U = 0.02;
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

const double T_ = 300;
const int N_step = T_/dt_;

class Cumu {
public:
    double9_t *f, *fpost, *fprev, *c, *cpost;
    int imax, jmax;
    double omega;

    Cumu(int imax, int jmax, double omega) : imax(imax), jmax(jmax), omega(omega) {
        f = new double9_t[imax*jmax];
        fpost = new double9_t[imax*jmax];
        fprev = new double9_t[imax*jmax];
        c = new double9_t[imax*jmax];
        cpost = new double9_t[imax*jmax];

        #pragma acc enter data \
        copyin(this[0:1], f[:imax*jmax], fpost[:imax*jmax], fprev[:imax*jmax], c[:imax*jmax], cpost[:imax*jmax])
    }

    ~Cumu() {
        delete[] f;
        delete[] fpost;
        delete[] fprev;
        delete[] c;
        delete[] cpost;

        #pragma acc exit data \
        delete(f[:imax*jmax], fpost[:imax*jmax], fprev[:imax*jmax], c[:imax*jmax], cpost[:imax*jmax], this[0:1])
    }

    void print_info() {
        printf("CUMULANT LBM\n");
        printf("\tdomain size = (%d %d)\n", imax, jmax);
        printf("\tghost cell = %d\n", Ghost_cell);
        printf("\trelaxation rate = %lf\n", omega);
    }
};

double9_t pdf_to_raw(const double9_t &f) {
    return double9_t{{
        f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8],
        -f[0] + f[2] - f[3] + f[5] - f[6] + f[8],
        f[0] + f[2] + f[3] + f[5] + f[6] + f[8],
        -f[0] - f[1] - f[2] + f[6] + f[7] + f[8],
        f[0] - f[2] - f[6] + f[8],
        -f[0] - f[2] + f[6] + f[8],
        f[0] + f[1] + f[2] + f[6] + f[7] + f[8],
        -f[0] + f[2] - f[6] + f[8],
        f[0] + f[2] + f[6] + f[8]
    }};
}

double9_t raw_to_pdf(const double9_t &m) {
    return double9_t{{
        0.25*(m[4] - m[5] - m[7] + m[8]),
        -0.5*(m[3] - m[5] - m[6] + m[8]),
        -0.25*(m[4] + m[5] - m[7] - m[8]),
        -0.5*(m[1] - m[2] - m[7] + m[8]),
        1.0*(m[0] - m[2] - m[6] + m[8]),
        0.5*(m[1] + m[2] - m[7] - m[8]),
        -0.25*(m[4] - m[5] + m[7] - m[8]),
        0.5*(m[3] - m[5] + m[6] - m[8]),
        0.25*(m[4] + m[5] + m[7] + m[8])
    }};
}

double9_t raw_to_center(const double9_t &m, const double u, const double v) {
    const double uu = u*u;
    const double vv = v*v;
    return double9_t{{
        m[0],
        -m[0]*v + m[1],
        m[0]*vv - 2*m[1]*v + m[2],
        -m[0]*u + m[3],
       -m[3]*v + m[4] + u*(m[0]*v - m[1]),
       m[3]*vv - 2*m[4]*v + m[5] - u*(m[0]*vv - 2*m[1]*v + m[2]),
       m[0]*uu - 2*m[3]*u + m[6],
       -m[6]*v + m[7] + uu*(-m[0]*v + m[1]) + 2*u*(m[3]*v - m[4]),
       m[6]*vv - 2*m[7]*v + m[8] + uu*(m[0]*vv - 2*m[1]*v + m[2]) - 2*u*(m[3]*vv - 2*m[4]*v + m[5])
    }};
}

double9_t center_to_raw(const double9_t &k, const double u, const double v) {
    const double uu = u*u;
    const double vv = v*v;
    return double9_t{{
        k[0],
        k[0]*v + k[1],
        k[0]*vv + 2*k[1]*v + k[2],
        k[0]*u + k[3],
        k[3]*v + k[4] + u*(k[0]*v + k[1]),
        k[3]*vv + 2*k[4]*v + k[5] + u*(k[0]*vv + 2*k[1]*v + k[2]),
        k[0]*uu + 2*k[3]*u + k[6],
        k[6]*v + k[7] + uu*(k[0]*v + k[1]) + 2*u*(k[3]*v + k[4]),
        k[6]*vv + 2*k[7]*v + k[8] + uu*(k[0]*vv + 2*k[1]*v + k[2]) + 2*u*(k[3]*vv + 2*k[4]*v + k[5])
    }};
}

void get_statistics(
    const double9_t &f,
    double &fluxx,
    double &fluxy,
    double &density
) {
    fluxx = 0;
    fluxy = 0;
    density = 0;
    for (int fid = 0; fid < 9; fid ++) {
        fluxx += Ve[fid][0]*f[fid];
        fluxy += Ve[fid][1]*f[fid];
        density += f[fid];
    }
}

void compute_cumulants(
    const double9_t *f,
    double9_t *c,
    int imax, 
    int jmax
) {
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], c[:imax*jmax], Ve, Cmk) \
    firstprivate(imax, jmax)
    for (int lat = 0; lat < imax*jmax; lat ++) {
        const double9_t &fc = f[lat];
        double fluxx, fluxy, density;
        get_statistics(fc, fluxx, fluxy, density);
        const double u = fluxx/density;
        const double v = fluxy/density;
        double9_t mc = pdf_to_raw(fc);
        double9_t kc = raw_to_center(mc, u, v);
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

double9_t get_equilibrium(double density, double u, double v) {
    return double9_t{{
        /* 00       01          02          10     11 12      20      21 22 */
        density, density*v, csq*density, density*u, 0, 0, csq*density, 0, 0
    }};
}

void relax_cumulants(
    const double9_t *c,
    double9_t *cpost,
    double omega,
    int imax,
    int jmax
) {
    #pragma acc parallel loop independent \
    present(c[:imax*jmax], cpost[:imax*jmax]) \
    firstprivate(omega, imax, jmax)
    for (int lat = 0; lat < imax*jmax; lat ++) {
        const double9_t &cc = c[lat];
        double9_t &ccpost = cpost[lat];
        ccpost[0] = cc[0];
        ccpost[1] = cc[1];
        ccpost[3] = cc[3];
        ccpost[2] = 0.5*(1 - omega)*(cc[2] - cc[6]) + cc[0]*csq;
        ccpost[4] = (1 - omega)*cc[4];
        ccpost[6] = 0.5*(1 - omega)*(cc[6] - cc[2]) + cc[0]*csq;
        ccpost[5] = 0;
        ccpost[7] = 0;
        ccpost[8] = 0;
    }
}

void compute_post_pdfs(
    const double9_t *cpost,
    double9_t *fpost,
    int imax,
    int jmax
) {
    #pragma acc parallel loop independent \
    present(cpost[:imax*jmax], fpost[:imax*jmax])\
    firstprivate(imax, jmax)
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

        double9_t mc = center_to_raw(kc, u, v);
        fpost[lat] = raw_to_pdf(mc);
    }
}

void advect(
    const double9_t *fpost,
    double9_t *f,
    int imax,
    int jmax
) {
    #pragma acc parallel loop independent collapse(3) \
    present(fpost[:imax*jmax], f[:imax*jmax], Ve) \
    firstprivate(imax, jmax)
    for (int i = 1; i < imax - 1; i ++) {
    for (int j = 1; j < jmax - 1; j ++) {
    for (int q = 0; q < 9; q ++) {
        int src = (i - Ve[q][0])*jmax + (j - Ve[q][1]);
        int dst = i*jmax + j;
        f[dst][q] = fpost[src][q];
    }}}
}

template<typename T>
void cpy_array(T *dst, T *src, int cnt) {
    #pragma acc parallel loop independent \
    present(dst[:cnt], src[:cnt]) \
    firstprivate(cnt)
    for (int i = 0; i < cnt; i ++) {
        dst[i] = src[i];
    }
}

void apply_fbc(
    double9_t *f,
    double9_t *fpost,
    double9_t *fprev,
    double u_inflow,
    int imax,
    int jmax
) {
    /* bottom slip */
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], fpost[:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int i = 1; i < imax - 1; i ++) {
        const int j = 1;
        const int lat = i*jmax + j;
        const int flist[]{2,5,8};
        for (int fid : flist) {
            const int ref = fid - 2;
            const int src = (i - Ve[ref][0])*jmax + j;
            f[lat][fid] = fpost[src][ref];
        }
    }
    /* top slip */
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], fpost[:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int i = 1; i < imax - 1; i ++) {
        const int j = jmax - 2;
        const int lat = i*jmax + j;
        const int flist[]{0,3,6};
        for (int fid : flist) {
            const int ref = fid + 2;
            const int src = (i - Ve[ref][0])*jmax + j;
            f[lat][fid] = fpost[src][ref];
        }
    }
    /* right convective outflow */
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], fprev[:imax*jmax], Ve, Wgt) \
    firstprivate(imax, jmax, u_inflow)
    for (int j = 1; j < jmax - 1; j ++) {
        const int i = imax - 2;
        const int lat = i*jmax + j;
        const int li1 = (i - 1)*jmax + j;
        const int li2 = (i - 2)*jmax + j;
        double u0, u1, u2, v0, v1, v2, rho0, rho1, rho2;
        get_statistics(fprev[lat], u0, v0, rho0);
        u0 /= rho0;
        v0 /= rho0;
        get_statistics(fprev[li1], u1, v1, rho1);
        u1 /= rho1;
        v1 /= rho1;
        get_statistics(fprev[li2], u2, v2, rho2);
        u2 /= rho2;
        v2 /= rho2;
        double2_t gradient{{
            0.5*(3*u0 - 4*u1 + u2),
            0.5*(3*v0 - 4*v1 + v2)
        }};
        const int flist[]{0,1,2};
        const double cs = sqrt(csq);
        for (int fid : flist) {
            f[lat][fid] = fprev[lat][fid] - 3*Wgt[fid]*u_inflow*rho0*(gradient[0]*Ve[fid][0] + gradient[1]*Ve[fid][1]);
        }
    }
    /* left inflow */
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], fpost[:imax*jmax], Ve, Wgt) \
    firstprivate(imax, jmax, u_inflow)
    for (int j = 1; j < jmax - 1; j ++) {
        const int i = 1;
        const int lat = i*jmax + j;
        const int flist[]{6,7,8};
        for (int fid : flist) {
            const int lnk = 8 - fid;
            f[lat][fid] = fpost[lat][lnk] - Ve[lnk][0]*2*u_inflow*Wgt[lnk]*csqi;
        }
    }

    const int ia = (5 - 0.5*L_)*Cells_per_length + 1;
    const int ib = (5 + 0.5*L_)*Cells_per_length + 1;
    const int ja = 0.5*(Ly_ - L_)*Cells_per_length + 1;
    const int jb = 0.5*(Ly_ + L_)*Cells_per_length + 1;
    /* square left face */
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], fpost[:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int j = ja; j < jb; j ++) {
        const int lat = (ia - 1)*jmax + j;
        const int flist[]{0,1,2};
        for (int fid : flist) {
            f[lat][fid] = fpost[lat][8 - fid];
        }
    }
    /* square right face */
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], fpost[:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int j = ja; j < jb; j ++) {
        const int lat = ib*jmax + j;
        const int flist[]{6,7,8};
        for (int fid : flist) {
            f[lat][fid] = fpost[lat][8 - fid];
        }
    }
    /* square bottom face */
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], fpost[:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int i = ia; i < ib; i ++) {
        const int lat = i*jmax + (ja - 1);
        const int flist[]{0,3,6};
        for (int fid : flist) {
            f[lat][fid] = fpost[lat][8 - fid];
        }
    }
    /* square top face */
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], fpost[:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int i = ia; i < ib; i ++) {
        const int lat = i*jmax + jb;
        const int flist[]{2,5,8};
        for (int fid : flist) {
            f[lat][fid] = fpost[lat][8 - fid];
        }
    }

    #pragma acc serial \
    present(f[:imax*jmax], fpost[:imax*jmax]) \
    firstprivate(imax, jmax)
    {
        f[(ia-1)*jmax + (ja-1)][0] = fpost[(ia-1)*jmax + (ja-1)][8];
        f[(ia-1)*jmax + (jb  )][2] = fpost[(ia-1)*jmax + (jb  )][6];
        f[(ib  )*jmax + (ja-1)][6] = fpost[(ib  )*jmax + (ja-1)][2];
        f[(ia  )*jmax + (jb  )][8] = fpost[(ia  )*jmax + (jb  )][0];
    }
}

Cumu *init(int imax, int jmax, double tau) {
    Cumu *cumu = new Cumu(imax, jmax, 1./tau);
    cumu->print_info();

    #pragma acc parallel loop independent \
    present(cumu[0:1], cumu->c[0:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int lat = 0; lat < imax*jmax; lat ++) {
        cumu->c[lat] = get_equilibrium(1, 0, 0);
    }

    compute_post_pdfs(cumu->c, cumu->f, cumu->imax, cumu->jmax);
    cpy_array(cumu->fpost, cumu->f, cumu->imax*cumu->jmax);
    cpy_array(cumu->fprev, cumu->f, cumu->imax*cumu->jmax);
    apply_fbc(cumu->f, cumu->fpost, cumu->fprev, U, cumu->imax, cumu->jmax);
    return cumu;
}

void finalize(Cumu *cumu) {
    delete cumu;
}

void main_loop(Cumu *cumu) {
    cpy_array(cumu->fprev, cumu->f, cumu->imax*cumu->jmax);
    compute_cumulants(cumu->f, cumu->c, cumu->imax, cumu->jmax);
    relax_cumulants(cumu->c, cumu->cpost, cumu->omega, cumu->imax, cumu->jmax);
    compute_post_pdfs(cumu->cpost, cumu->fpost, cumu->imax, cumu->jmax);
    advect(cumu->fpost, cumu->f, cumu->imax, cumu->jmax);
    apply_fbc(cumu->f, cumu->fpost, cumu->fprev, U, cumu->imax, cumu->jmax);
}

void output(Cumu *cumu) {
    #pragma acc update \
    self(cumu->f[:cumu->imax*cumu->jmax])

    FILE *file = fopen("data/cumu.csv", "w");
    fprintf(file, "x,y,z,u,v,w,rho\n");
    for (int j = 1; j < cumu->jmax - 1; j ++) {
    for (int i = 1; i < cumu->imax - 1; i ++) {
        const int lat = i*cumu->jmax + j;
        double density, fluxx, fluxy;
        get_statistics(cumu->f[lat], fluxx, fluxy, density);
        double u = fluxx/density, v = fluxy/density;
        fprintf(
            file,
            "%e,%e,%e,%lf,%lf,%lf,%lf\n",
            (i-1+0.5)*dx_, (j-1+0.5)*dx_, 0.0,
            u*Cu_, v*Cu_, 0.0,
            density
        );
    }}
    fclose(file);
}

void output(Cumu *cumu, int n) {
    #pragma acc update \
    self(cumu->f[:cumu->imax*cumu->jmax])

    FILE *file = fopen(("data/cumu.csv." + std::to_string(n)).c_str(), "w");
    fprintf(file, "x,y,z,u,v,w,rho\n");
    for (int j = 1; j < cumu->jmax - 1; j ++) {
    for (int i = 1; i < cumu->imax - 1; i ++) {
        const int lat = i*cumu->jmax + j;
        double density, fluxx, fluxy;
        get_statistics(cumu->f[lat], fluxx, fluxy, density);
        double u = fluxx/density, v = fluxy/density;
        fprintf(
            file,
            "%e,%e,%e,%lf,%lf,%lf,%lf\n",
            (i-1+0.5)*dx_, (j-1+0.5)*dx_, 0.0,
            u*Cu_, v*Cu_, 0.0,
            density
        );
    }}
    fclose(file);
}

int main() {
    Cumu *cumu = init(Lx_*Cells_per_length + 2*Ghost_cell, Ly_*Cells_per_length + 2*Ghost_cell, tau);
    auto begin = std::chrono::high_resolution_clock::now();
    for (int step = 1; step <= N_step; step ++) {
        main_loop(cumu);
        printf("\r%d/%d", step, N_step);
        fflush(stdout);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapse = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    printf("\n%d\n", elapse);
    output(cumu, 0);
    for (int step = 1; step <= int(100/dt_); step ++) {
        main_loop(cumu);
        printf("\r%d/%d", step, int(100/dt_));
        fflush(stdout);
        if (step % int(1/dt_) == 0) {
            output(cumu, step / int(1/dt_));
        }
    }
    printf("\n");
    finalize(cumu);
}
