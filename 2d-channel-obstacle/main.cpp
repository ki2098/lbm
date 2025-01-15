#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <chrono>
#include <iostream>
#include <fstream>

template<typename T, int N>
struct vector_t {
    T m[N];

    T &operator[](int i) {
        return m[i];
    }
    
    const T &operator[](int i) const {
        return m[i];
    }
    std::string to_str() {
        std::string str = "(";
        if (N == 0) {
            return str + ")";
        }
        str += std::to_string(m[0]);
        for (int i = 1; i < N; i ++) {
            str += ", " + std::to_string(m[i]);
        }
        return str + ")";
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
const double Lx_ = 50*L_;
const double Ly_ = 25*L_;
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
const double Cpressure_ = 1.*Cu_*Cu_;

const double T_pre = 500;
const int N_pre = T_pre/dt_;
const double T_post = 500;
const int N_post = T_post/dt_;

const double center_x_ = 10*L_;
const double center_y_ = Ly_/2 + 2*L_;

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
        printf("\tRe = %lf\n", Re);
    }
};

void get_statistics(
    const double9_t &f,
    double &phix,
    double &phiy,
    double &rho
) {
    phix = 0;
    phiy = 0;
    rho = 0;
    for (int fid = 0; fid < 9; fid ++) {
        phix += Ve[fid][0]*f[fid];
        phiy += Ve[fid][1]*f[fid];
        rho += f[fid];
    }
}


double9_t get_equilibrium_cumulant(double rho, double u, double v) {
    return double9_t{{
        /* 00       01          02          10     11 12      20      21 22 */
        rho, rho*v, csq*rho, rho*u, 0, 0, csq*rho, 0, 0
    }};
}

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

double9_t raw_to_central(const double9_t &m, const double u, const double v) {
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

double9_t central_to_raw(const double9_t &k, const double u, const double v) {
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

double9_t central_to_cumulant(const double9_t &k, const double u, const double v) {
    double9_t c;
    const double rho = k[0];
    c[0] = rho;
    c[1] = v*rho;
    c[3] = u*rho;
    c[2] = k[2];
    c[4] = k[4];
    c[6] = k[6];
    c[5] = k[5];
    c[7] = k[7];
    c[8] = k[8] - (2*k[4]*k[4] + k[6]*k[2])/rho;
    return c;
}


double9_t cumulant_to_central(const double9_t &c) {
    double9_t k;
    const double rho = c[0];
    k[0] = rho;
    k[1] = 0;
    k[3] = 0;
    k[2] = c[2];
    k[4] = c[4];
    k[6] = c[6];
    k[5] = c[5];
    k[7] = c[7];
    k[8] = c[8] + (2*k[4]*k[4] + k[6]*k[2])/rho;
    return k;
}

double9_t pdf_to_cumulant(const double9_t &f) {
    double phix, phiy, rho;
    get_statistics(f, phix, phiy, rho);
    const double u = phix/rho;
    const double v = phiy/rho;
    return central_to_cumulant(raw_to_central(pdf_to_raw(f), u, v), u, v);
}

double9_t cumulant_to_pdf(const double9_t &c) {
    double rho = c[0];
    double v = c[1]/rho;
    double u = c[3]/rho;
    return raw_to_pdf(central_to_raw(cumulant_to_central(c), u, v));
}

void compute_cumulant(
    const double9_t *f,
    double9_t *c,
    int imax, 
    int jmax
) {
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], c[:imax*jmax], Ve, Cmk) \
    firstprivate(imax, jmax)
    for (int lat = 0; lat < imax*jmax; lat ++) {
        c[lat] = pdf_to_cumulant(f[lat]);
    }
}

double9_t relax_cumulant(const double9_t &c, const double omega) {
    double9_t cpost;
    const double rho = c[0];
    cpost[0] = c[0];
    cpost[1] = c[1];
    cpost[3] = c[3];
    cpost[2] = 0.5*(1 - omega)*(c[2] - c[6]) + rho*csq;
    cpost[4] = (1 - omega)*c[4];
    cpost[6] = 0.5*(1 - omega)*(c[6] - c[2]) + rho*csq;
    cpost[5] = 0;
    cpost[7] = 0;
    cpost[8] = 0;
    return cpost;
}

void collide(
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
        cpost[lat] = relax_cumulant(c[lat], omega);
    }
}

void compute_post_collision_pdfs(
    const double9_t *cpost,
    double9_t *fpost,
    int imax,
    int jmax
) {
    #pragma acc parallel loop independent \
    present(cpost[:imax*jmax], fpost[:imax*jmax])\
    firstprivate(imax, jmax)
    for (int lat = 0; lat < imax*jmax; lat ++) {
        fpost[lat] = cumulant_to_pdf(cpost[lat]);
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
        double ru0, ru1, ru2, rv0, rv1, rv2, dummy;
        get_statistics(fprev[lat], ru0, rv0, dummy);
        get_statistics(fprev[li1], ru1, rv1, dummy);
        get_statistics(fprev[li2], ru2, rv2, dummy);
        double2_t gradient{{
            (3*ru0 - 4*ru1 + ru2),
            (3*rv0 - 4*rv1 + rv2)
        }};
        const int flist[]{0,1,2};
        const double cs = sqrt(csq);
        for (int fid : flist) {
            f[lat][fid] = fprev[lat][fid] - 3*Wgt[fid]*0.5*u_inflow*(gradient[0]*Ve[fid][0] + gradient[1]*Ve[fid][1]);
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
        // const int i = 0;
        // const int lat = i*jmax + j;
        // f[lat] = cumulant_to_pdf(get_equilibrium_cumulant(1., u_inflow, 0.));
    }

    const int ia = (center_x_ - 0.5*L_)*Cells_per_length + 1;
    const int ib = (center_x_ + 0.5*L_)*Cells_per_length + 1;
    const int ja = (center_y_ - 0.5*L_)*Cells_per_length + 1;
    const int jb = (center_y_ + 0.5*L_)*Cells_per_length + 1;
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
        f[(ia - 1)*jmax + (ja - 1)][0] = fpost[(ia - 1)*jmax + (ja - 1)][8];
        f[(ia - 1)*jmax + (jb    )][2] = fpost[(ia - 1)*jmax + (jb    )][6];
        f[(ib    )*jmax + (ja - 1)][6] = fpost[(ib    )*jmax + (ja - 1)][2];
        f[(ib    )*jmax + (jb    )][8] = fpost[(ib    )*jmax + (jb    )][0];
    }
}

class CdMonitor {
public:
    int imax, jmax;
    int ia, ib, ja, jb;
    bool file_o;
    std::ofstream ofs;

    CdMonitor(int ia, int ib, int ja, int jb, int imax, int jmax, bool file_o):
    ia(ia), ib(ib), ja(ja), jb(jb), imax(imax), jmax(jmax), file_o(file_o) {
        if (file_o) {
            ofs.open("data/cd_history.csv");
            ofs << "t,cd\n";
        }
        #pragma acc enter data \
        copyin(this[0:1])
    }

    double get_cd(double9_t *f, double t) {
        double drag = 0;
        /* pressure drag on left face */
        #pragma acc parallel loop independent \
        reduction(+:drag) \
        present(f[:imax*jmax], this[0:1]) 
        for (int j = ja; j < jb; j ++) {
            const int i = ia - 1;
            const int lat = i*jmax + j;
            double phix, phiy, rho;
            get_statistics(f[lat], phix, phiy, rho);
            drag += (rho - 1)*csq*dx;
        }

        /* pressure drag on right face */
        #pragma acc parallel loop independent \
        reduction(+:drag) \
        present(f[:imax*jmax], this[0:1]) 
        for (int j = ja; j < jb; j ++) {
            const int i = ib;
            const int lat = i*jmax + j;
            double phix, phiy, rho;
            get_statistics(f[lat], phix, phiy, rho);
            drag += (- (rho - 1)*csq*dx);
        }

        /* shear drag on upper face */
        #pragma acc parallel loop independent \
        reduction(+:drag) \
        present(f[:imax*jmax], this[0:1]) 
        for (int i = ia; i < ib; i ++) {
            const int j = jb;
            const int lat = i*jmax + j;
            double phix, phiy, rho;
            get_statistics(f[lat], phix, phiy, rho);
            drag += 2*nu*phix*dx;
        }

        /* shear drag on lower face */
        #pragma acc parallel loop independent \
        reduction(+:drag) \
        present(f[:imax*jmax], this[0:1]) 
        for (int i = ia; i < ib; i ++) {
            const int j = ja - 1;
            const int lat = i*jmax + j;
            double phix, phiy, rho;
            get_statistics(f[lat], phix, phiy, rho);
            drag += 2*nu*phix*dx;
        }

        double cd = 2*drag/(1.*U*U*(jb - ja)*dx);
        if (file_o) {
            ofs << t << "," << cd << std::endl;
        }
        return cd;
    }

    ~CdMonitor() {
        #pragma acc exit data \
        delete(this[0:1])
        if (file_o) {
            ofs.close();
        }
    }
};

Cumu *init(int imax, int jmax, double tau) {
    Cumu *cumu = new Cumu(imax, jmax, 1./tau);
    cumu->print_info();

     #pragma acc parallel loop independent \
    present(cumu[0:1], cumu->f[0:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int lat = 0; lat < imax*jmax; lat ++) {
        cumu->f[lat] = cumulant_to_pdf(get_equilibrium_cumulant(1, U, 0));
    }
    return cumu;
}

void finalize(Cumu *cumu) {
    delete cumu;
}

void main_loop(Cumu *cumu) {
    cpy_array(cumu->fprev, cumu->f, cumu->imax*cumu->jmax);
    compute_cumulant(cumu->f, cumu->c, cumu->imax, cumu->jmax);
    collide(cumu->c, cumu->cpost, cumu->omega, cumu->imax, cumu->jmax);
    compute_post_collision_pdfs(cumu->cpost, cumu->fpost, cumu->imax, cumu->jmax);
    advect(cumu->fpost, cumu->f, cumu->imax, cumu->jmax);
    apply_fbc(cumu->f, cumu->fpost, cumu->fprev, U, cumu->imax, cumu->jmax);
}

void output(Cumu *cumu) {
    #pragma acc update \
    self(cumu->f[:cumu->imax*cumu->jmax])

    FILE *file = fopen("data/cumu.csv", "w");
    fprintf(file, "x,y,z,u,v,w,p\n");
    for (int j = 1; j < cumu->jmax - 1; j ++) {
    for (int i = 1; i < cumu->imax - 1; i ++) {
        const int lat = i*cumu->jmax + j;
        double rho, phix, phiy;
        get_statistics(cumu->f[lat], phix, phiy, rho);
        double u = phix/rho, v = phiy/rho, p = (rho - 1)*csq*Cpressure_;
        fprintf(
            file,
            "%e,%e,%e,%lf,%lf,%lf,%lf\n",
            (i - 1 + 0.5)*dx_, (j - 1 + 0.5)*dx_, 0.0,
            u*Cu_, v*Cu_, 0.0,
            p
        );
    }}
    fclose(file);
}

void output(Cumu *cumu, int n) {
    #pragma acc update \
    self(cumu->f[:cumu->imax*cumu->jmax])

    FILE *file = fopen(("data/cumu.csv." + std::to_string(n)).c_str(), "w");
    fprintf(file, "x,y,z,u,v,w,p\n");
    for (int j = 1; j < cumu->jmax - 1; j ++) {
    for (int i = 1; i < cumu->imax - 1; i ++) {
        const int lat = i*cumu->jmax + j;
        double rho, phix, phiy;
        get_statistics(cumu->f[lat], phix, phiy, rho);
        double u = phix/rho, v = phiy/rho, p = (rho - 1)*csq*Cpressure_;
        fprintf(
            file,
            "%e,%e,%e,%lf,%lf,%lf,%lf\n",
            (i - 1 + 0.5)*dx_, (j - 1 + 0.5)*dx_, 0.0,
            u*Cu_, v*Cu_, 0.0,
            p
        );
    }}
    fclose(file);
}

class Probe {
public:
    int i, j, imax, jmax;
    std::ofstream ofs;

    Probe(int i, int j, int imax, int jmax) : i(i), j(j), imax(imax), jmax(jmax) {
        ofs.open("data/probe.csv");
        ofs << "t,value\n";
    }

    void write(const double9_t *f, double t) {
        int lat = i*jmax + j;
        #pragma acc update \
        self(f[lat:1])

        double rv, ru, r;
        get_statistics(f[lat], ru, rv, r);
        ofs << t << "," << ru/r << std::endl;
    }

    ~Probe() {
        ofs.close();
    }

};

int main() {  
    Cumu *cumu = init(Lx_*Cells_per_length + 2*Ghost_cell, Ly_*Cells_per_length + 2*Ghost_cell, tau);
    Probe probe(
        (center_x_ + L_)*Cells_per_length + 1,
        (center_y_ + 0.5*L_)*Cells_per_length + 1,
        Lx_*Cells_per_length + 2*Ghost_cell,
        Ly_*Cells_per_length + 2*Ghost_cell
    );
    CdMonitor cd_monitor(
        (center_x_ - 0.5*L_)*Cells_per_length + 1,
        (center_x_ + 0.5*L_)*Cells_per_length + 1,
        (center_y_ - 0.5*L_)*Cells_per_length + 1,
        (center_y_ + 0.5*L_)*Cells_per_length + 1,
        Lx_*Cells_per_length + 2*Ghost_cell,
        Ly_*Cells_per_length + 2*Ghost_cell,
        true
    );
    auto begin = std::chrono::high_resolution_clock::now();
    for (int step = 1; step <= N_pre; step ++) {
        main_loop(cumu);
        printf("\r%d/%d", step, N_pre);
        fflush(stdout);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapse = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    printf("\n%d\n", elapse);

    for (int step = 1; step <= N_post; step ++) {
        main_loop(cumu);
        printf("\r%d/%d", step, N_post);
        fflush(stdout);
        if (step % int(1/dt_) == 0) {
            // output(cumu, step / int(1/dt_));
        }
        probe.write(cumu->f, (step + N_pre)*dt_);
        cd_monitor.get_cd(cumu->f, (step + N_pre)*dt_);
    }
    output(cumu);
    // drag_history.close();
    printf("\n");
    finalize(cumu);
}
