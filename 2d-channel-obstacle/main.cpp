#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>

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
const double Re = 1e5;
const double Lx_ = 20*L_;
const double Ly_ =  10*L_;
const int Cells_per_length = 20;
const int Ghost_cell = 1;
const double nu_ = L_*U_/Re;
const double dx_ = 1./Cells_per_length;
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
const double csqi = 3;
const double tau = nu*csqi + 0.5;
const double omega = 1/tau;
const double Re_lattice = U*dx/nu;

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
        double9_t kc{{}};
        for (int fid = 0; fid < 9; fid ++) {
            const int i = Ve[fid][0];
            const int j = Ve[fid][1];
            vector_t<double, 3> ca, cb;
            ca[0] = cb[0] = 1;
            for (int order = 1; order < 3; order ++) {
                ca[order] = ca[order - 1]*(i - u);
                cb[order] = cb[order - 1]*(j - v);
            }
            for (int kid = 0; kid < 9; kid ++) {
                const int a = Cmk[kid][0];
                const int b = Cmk[kid][1];
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
    #pragma acc parallel loop independent collapse(3) \
    present(fpost[:imax*jmax], f[:imax*jmax], Ve) \
    firstprivate(imax, jmax)
    for (int i = 1; i < imax - 1; i ++) {
    for (int j = 1; j < jmax - 1; j ++) {
    for (int q = 0; q < 9; q ++) {
        int dst = (i + Ve[q][0])*jmax + (j + Ve[q][1]);
        int src = i*jmax + j;
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
    /* right outflow */
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], fprev[:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int j = 1; j < jmax - 1; j ++) {
        const int i = imax - 2;
        const int lat = i*jmax + j;
        const int lin = (i - 1)*jmax + j;
        const int flist[]{0,1,2};
        const double cs = sqrt(csq);
        for (int fid : flist) {
            f[lat][fid] = cs*fprev[lin][fid] + (1 - cs)*fprev[lat][fid];
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
}

Cumu *init(int imax, int jmax, double tau) {
    Cumu *cumu = new Cumu(imax, jmax, 1./tau);
    cumu->print_info();

    #pragma acc parallel loop independent \
    present(cumu[0:1], cumu->c[0:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int lat = 0; lat < imax*jmax; lat ++) {
        cumu->c[lat] = get_equilibrium(1, U, 0);
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
    for (int step = 1; step <= N_step; step ++) {
        main_loop(cumu);
        printf("\r%d/%d", step, N_step);
        fflush(stdout);
    }
    printf("\n");
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
    // output(cumu);
    finalize(cumu);
}
