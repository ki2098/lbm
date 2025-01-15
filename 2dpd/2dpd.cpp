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
const double Re = 2.9e5;
const double Ox_ = -5*L_;
const double Oy_ = -7.5*L_;
const double Lx_ = 30*L_;
const double Ly_ = -2*Oy_;
const int Cells_per_length = 20;
const int Ghost_cell = 1;
const double nu_ = L_*U_/Re;
const double dx_ = 1./Cells_per_length;
const double dx = 1;
const double Cx_ = dx_/dx;
const double Ox = Ox_/Cx_;
const double Oy = Oy_/Cx_;
const double U = 0.025;
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

const double CDD_ = 13;
const double Turbine_xy_[2][2] = {
    {0., 0.},
    {5., 1.}
};
const double Diameter_ = 1.5;
const double Thick_ = 0.03;

const double Ccdd_ = 1./Cx_;
const double CDD = CDD_/Ccdd_;
const double Turbine_xy[2][2] = {
    {Turbine_xy_[0][0]/Cx_, Turbine_xy_[0][1]/Cx_},
    {Turbine_xy_[1][0]/Cx_, Turbine_xy_[1][1]/Cx_}
};
#pragma acc declare \
copyin(Turbine_xy)
const double Diameter = Diameter_/Cx_;
const double Thick = Thick_/Cx_;

const double T_ = 300;
const int Nt = T_/dt_;

template<typename T>
T sq(T a) {
    return a*a;
}

double intersection_length(double x0, double x1, double x2, double x3) {
    if (x1 > x2 && x0 < x3) {
        return fmin(x3-x2, fmin(x1-x0, fmin(x1-x2, x3-x0)));
    } else {
        return 0;
    }
}

class PorousDisk {
public:
    double *dfunc;
    int imax, jmax;

    PorousDisk(int imax, int jmax) : imax(imax), jmax(jmax) {
        dfunc = new double[imax*jmax]{};
        for (int i = 0; i < imax; i ++) {
        for (int j = 0; j < jmax; j ++) {
            double x0 = Ox + dx*(i - 1);
            double x1 = Ox + dx*(i);
            double y0 = Oy + dx*(j - 1);
            double y1 = Oy + dx*(j);
            for (int t = 1; t < 2; t ++) {
                double x2 = Turbine_xy[t][0] - 0.5*Thick;
                double x3 = Turbine_xy[t][0] + 0.5*Thick;
                double y2 = Turbine_xy[t][1] - 0.5*Diameter;
                double y3 = Turbine_xy[t][1] + 0.5*Diameter;
                double cover = intersection_length(x0, x1, x2, x3)*intersection_length(y0, y1, y2, y3);
                if (cover > 0) {
                    double sigma = Diameter/6;
                    double d = fabs(0.5*(y0 + y1) - Turbine_xy[t][1]);
                    const int lat = i*jmax + j;
                    dfunc[lat] = CDD*exp(-0.5*sq(d/sigma))*cover/(dx*dx);
                    break;
                }
            }
        }}

        #pragma acc enter data \
        copyin(this[0:1], dfunc[:imax*jmax])
    }

    ~PorousDisk() {
        delete[] dfunc;

        #pragma acc exit data \
        delete(dfunc[:imax*jmax], this[0:1])
    }

    void print_info() {
        printf("PD INFO\n");
        printf("\tGlobal size = (%d %d)\n", imax, jmax);
        printf("\tCDD = %lf\n", CDD_);
    }

    void output_dfunc(std::string path) {
        FILE *file = fopen(path.c_str(), "w");
        fprintf(file, "x,y,z,dfunc\n");
        for (int j = 1; j < jmax - 1; j ++) {
        for (int i = 1; i < imax - 1; i ++) {
            fprintf(
                file,
                "%e,%e,%e,%e\n",
                (i - 0.5)*dx_ + Ox_, (j - 0.5)*dx_ + Oy_, 0.0,
                dfunc[i*jmax + j]*Ccdd_
            );
        }}
        fclose(file);
    }
};

class Cumu {
public:
    double9_t *f, *fpost, *fprev, *c, *cpost;
    double2_t *shift;
    double3_t *mac; /* [ρ u v] */

    int imax, jmax;
    double omega;

    Cumu(int imax, int jmax, double omega) : imax(imax), jmax(jmax), omega(omega) {
        f = new double9_t[imax*jmax];
        fpost = new double9_t[imax*jmax];
        fprev = new double9_t[imax*jmax];
        c = new double9_t[imax*jmax];
        cpost = new double9_t[imax*jmax];
        shift = new double2_t[imax*jmax];
        mac = new double3_t[imax*jmax];

        #pragma acc enter data \
        copyin(this[0:1], f[:imax*jmax], fpost[:imax*jmax], fprev[:imax*jmax], c[:imax*jmax], cpost[:imax*jmax], shift[:imax*jmax], mac[:imax*jmax])
    }

    ~Cumu() {
        delete[] f;
        delete[] fpost;
        delete[] fprev;
        delete[] c;
        delete[] cpost;
        delete[] shift;
        delete[] mac;

        #pragma acc exit data \
        delete(f[:imax*jmax], fpost[:imax*jmax], fprev[:imax*jmax], c[:imax*jmax], cpost[:imax*jmax], shift[:imax*jmax], mac[:imax*jmax], this[0:1])
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
    for (int q = 0; q < 9; q ++) {
        density += f[q];
        u += Ve[q][0]*f[q];
        v += Ve[q][1]*f[q];
    }
    u = (u + shift[0])/density; /* u = (m10 + φx/2)/m00 */
    v = (v + shift[1])/density; /* v = (m01 + φy/2)/m00 */
    return double3_t{{density, u, v}};
}

double9_t get_equilibrium_cumulant(const double3_t &mac) {
    double density = mac[0];
    double u       = mac[1];
    double v       = mac[2];
    return double9_t{{
        /* 00       01          02          10     11 12      20      21 22 */
        density, density*v, csq*density, density*u, 0, 0, csq*density, 0, 0
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

double9_t raw_to_central(const double9_t &m, const double3_t &mac) {
    const double u  = mac[1];
    const double v  = mac[2];
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

double9_t central_to_raw(const double9_t &k, const double3_t &mac) {
    const double u  = mac[1];
    const double v  = mac[2];
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

double9_t central_to_cumulant(const double9_t &k, const double3_t &mac) {
    const double density = k[0];
    const double u     = mac[1];
    const double v     = mac[2];
    double9_t c;
    c[0] = density;
    c[1] = v*density;
    c[3] = u*density;
    c[2] = k[2];
    c[4] = k[4];
    c[6] = k[6];
    c[5] = k[5];
    c[7] = k[7];
    c[8] = k[8] - (2*k[4]*k[4] + k[6]*k[2])/density;
    return c;
}

double9_t cumulant_to_central(const double9_t &c, const double2_t &shift) {
    const double density = c[0];
    const double sx  = shift[0];
    const double sy  = shift[1];
    double9_t k;
    k[0] = density;
    k[1] = sy; /* φy/2 */
    k[3] = sx; /* φx/2 */
    k[2] = c[2];
    k[4] = c[4];
    k[6] = c[6];
    k[5] = c[5];
    k[7] = c[7];
    k[8] = c[8] + (2*k[4]*k[4] + k[6]*k[2])/density;
    return k;
}

double9_t pdf_to_cumulant(const double9_t &f, const double3_t &mac) {
    return central_to_cumulant(raw_to_central(pdf_to_raw(f), mac), mac);
}

double9_t cumulant_to_pdf(const double9_t &c, const double3_t &mac, const double2_t &shift) {
    return raw_to_pdf(central_to_raw(cumulant_to_central(c, shift), mac));
}

void compute_cumulant(
    const double9_t *f,
    const double3_t *mac,
    double9_t *c,
    int imax, 
    int jmax
) {
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], mac[:imax*jmax], c[:imax*jmax], Ve, Cmk) \
    firstprivate(imax, jmax)
    for (int lat = 0; lat < imax*jmax; lat ++) {
        c[lat] = pdf_to_cumulant(f[lat], mac[lat]);
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

void compute_post_collision_pdf(
    const double9_t *cpost,
    const double3_t *mac,
    const double2_t *shift,
    double9_t *fpost,
    int imax,
    int jmax
) {
    #pragma acc parallel loop independent \
    present(cpost[:imax*jmax], mac[:imax*jmax], shift[:imax*jmax], fpost[:imax*jmax])\
    firstprivate(imax, jmax)
    for (int lat = 0; lat < imax*jmax; lat ++) {
        fpost[lat] = cumulant_to_pdf(cpost[lat], mac[lat], shift[lat]);
    }
}

void compute_macroscopic_field(
    const double9_t *f,
    const double2_t *shift,
    double3_t *mac,
    int imax,
    int jmax
) {
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], shift[:imax*jmax], mac[:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int lat = 0; lat < imax*jmax; lat ++) {
        mac[lat] = get_macroscopic_val(f[lat], shift[lat]);
    }
}

void compute_shift(
    const double3_t *mac,
    const double *dfunc,
    double2_t *shift,
    int imax,
    int jmax
) {
    #pragma acc parallel loop independent \
    present(mac[:imax*jmax], dfunc[:imax*jmax], shift[:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int lat = 0; lat < imax*jmax; lat ++) {
        const double u = mac[lat][1];
        const double v = mac[lat][2];
        const double Uabs = sqrt(u*u + v*v);
        shift[lat][0] = - 0.5*u*Uabs*dfunc[lat];
        shift[lat][1] = - 0.5*v*Uabs*dfunc[lat];
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
    const double9_t *fpost,
    const double9_t *fprev,
    const double3_t *mac,
    double9_t *f,
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
        const int qlist[]{2,5,8};
        for (int q : qlist) {
            const int ref = q - 2;
            const int src = (i - Ve[ref][0])*jmax + j;
            f[lat][q] = fpost[src][ref];
        }
    }
    /* top slip */
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], fpost[:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int i = 1; i < imax - 1; i ++) {
        const int j = jmax - 2;
        const int lat = i*jmax + j;
        const int qlist[]{0,3,6};
        for (int q : qlist) {
            const int ref = q + 2;
            const int src = (i - Ve[ref][0])*jmax + j;
            f[lat][q] = fpost[src][ref];
        }
    }
    // /* right convective outflow */
    // #pragma acc parallel loop independent \
    // present(f[:imax*jmax], fprev[:imax*jmax], mac[:imax*jmax], Ve, Wgt) \
    // firstprivate(imax, jmax, u_inflow)
    // for (int j = 1; j < jmax - 1; j ++) {
    //     // const int i = imax - 2;
    //     // const int lat = i*jmax + j;
    //     // const int li1 = (i - 1)*jmax + j;
    //     // const int li2 = (i - 2)*jmax + j;
    //     // double u0 = mac[lat][1];
    //     // double u1 = mac[li1][1];
    //     // double u2 = mac[li2][1];
    //     // double v0 = mac[lat][2];
    //     // double v1 = mac[li1][2];
    //     // double v2 = mac[li2][2];
    //     // double dudx = 0.5*(3*u0 - 4*u1 + u2);
    //     // double dvdx = 0.5*(3*v0 - 4*v1 + v2);
    //     // const int qlist[]{0,1,2};
    //     // const double cs = sqrt(csq);
    //     // for (int q : qlist) {
    //     //     f[lat][q] = fprev[lat][q] - 3*Wgt[q]*u0*(dudx*Ve[q][0] + dvdx*Ve[q][1]);
    //     // }
    //     const int i = imax - 2;
    //     const int lat = i*jmax + j;
    //     const int lin = (i - 1)*jmax + j;
    //     double u0 = mac[lat][1];
    //     double u1 = mac[lin][1];
    //     double v0 = mac[lat][2];
    //     double v1 = mac[lin][2];
    //     double dudx = u0 - u1;
    //     double dvdx = v0 - v1;
    //     const int qlist[]{0,1,2};
    //     const double cs = sqrt(csq);
    //     for (int q : qlist) {
    //         f[lat][q] = fprev[lat][q] - 3*Wgt[q]*u0*(dudx*Ve[q][0] + dvdx*Ve[q][1]);
    //     }
    // }
    /* right extrapolation outflow */
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], fprev[:imax*jmax], mac[:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int j = 1; j < jmax - 1; j ++) {
        const int i = imax - 2;
        const int lat = i*jmax + j;
        const int lin = (i - 1)*jmax + j;
        const int qlist[]{0, 1, 2};
        for (int q : qlist) {
            double mix = sqrt(csq) - mac[lat][1];
            f[lat][q] = mix*fprev[lin][q] + (1. - mix)*fprev[lat][q];
        }
    }
    /* left inflow */
    #pragma acc parallel loop independent \
    present(f[:imax*jmax], fpost[:imax*jmax], Ve, Wgt) \
    firstprivate(imax, jmax, u_inflow)
    for (int j = 1; j < jmax - 1; j ++) {
        const int i = 1;
        const int lat = i*jmax + j;
        const int qlist[]{6,7,8};
        for (int q : qlist) {
            const int lnk = 8 - q;
            f[lat][q] = fpost[lat][lnk] - Ve[lnk][0]*2*u_inflow*Wgt[lnk]*csqi;
        }
        // const int i = 0;
        // const int lat = i*jmax + j;
        // f[lat] = cumulant_to_pdf(get_equilibrium_cumulant(1., u_inflow, 0.));
    }
}

void init(int imax, int jmax, double tau, Cumu *&cumu, PorousDisk *&pd) {
    PorousDisk *_pd = new PorousDisk(imax, jmax);
    _pd->print_info();
    _pd->output_dfunc("data/dfunc.csv");

    Cumu *_cumu = new Cumu(imax, jmax, 1./tau);
    _cumu->print_info();
    #pragma acc parallel loop independent \
    present(_cumu[0:1], _cumu->mac[0:imax*jmax], _cumu->shift[:imax*jmax], _cumu->f[:imax*jmax]) \
    firstprivate(imax, jmax)
    for (int lat = 0; lat < imax*jmax; lat ++) {
        double3_t mac{{1., U, 0.}};
        double2_t shift{{0., 0.}};
        double9_t ceq = get_equilibrium_cumulant(mac);
        _cumu->mac[lat] = mac;
        _cumu->shift[lat] = shift;
        _cumu->f[lat] = cumulant_to_pdf(ceq, mac, shift);
    }
    compute_macroscopic_field(_cumu->f, _cumu->shift, _cumu->mac, imax, jmax);
    compute_shift(_cumu->mac, _pd->dfunc, _cumu->shift, _cumu->imax, _cumu->jmax);
    cumu = _cumu;
    pd = _pd;
}

void finalize(Cumu *cumu, PorousDisk *pd) {
    delete cumu;
    delete pd;
}

void output(Cumu *cumu, std::string path) {
    int imax = cumu->imax, jmax = cumu->jmax;
    #pragma acc update \
    self(cumu->mac[:imax*jmax])

    FILE *file = fopen(path.c_str(), "w");
    fprintf(file, "x,y,z,u,v,w,p\n");
    for (int j = 1; j < jmax - 1; j ++) {
    for (int i = 1; i < imax - 1; i ++) {
        const int lat = i*jmax + j;
        double density = cumu->mac[lat][0];
        double u = cumu->mac[lat][1];
        double v = cumu->mac[lat][2];
        double p = (density - 1.)*csq;
        fprintf(
            file,
            "%e,%e,%e,%lf,%lf,%lf,%lf\n",
            (i - 0.5)*dx_ + Ox_, (j - 0.5)*dx_ + Oy_, 0.,
            u*Cu_, v*Cu_, 0.,
            p*Cpressure_
        );
    }}
    fclose(file);
}

void main_loop(Cumu *cumu, PorousDisk *pd) {
    int imax = cumu->imax, jmax = cumu->jmax;
    cpy_array(cumu->fprev, cumu->f, imax*jmax);
    compute_cumulant(cumu->f, cumu->mac, cumu->c, imax, jmax);
    collide(cumu->c, cumu->cpost, cumu->omega, imax, jmax);
    compute_post_collision_pdf(cumu->cpost, cumu->mac, cumu->shift, cumu->fpost, imax, jmax);
    advect(cumu->fpost, cumu->f, imax, jmax);
    apply_fbc(cumu->fpost, cumu->fprev, cumu->mac, cumu->f, U, imax, jmax);
    compute_macroscopic_field(cumu->f, cumu->shift, cumu->mac, imax, jmax);
    compute_shift(cumu->mac, pd->dfunc, cumu->shift, imax, jmax);
}

int main() {
    const int imax = Lx_*Cells_per_length + 2*Ghost_cell;
    const int jmax = Ly_*Cells_per_length + 2*Ghost_cell;
    Cumu *cumu;
    PorousDisk *pd;
    init(imax, jmax, tau, cumu, pd);

    for (int step = 1; step <= Nt; step ++) {
        main_loop(cumu, pd);
        printf("\r%d/%d", step, Nt);
        fflush(stdout);

        if (step*dt_ >= T_ - 10) {
            if (step % int(0.5/dt_) == 0) {
                std::string path = "data/o.csv." + std::to_string(step / int(0.5/dt_));
                output(cumu, path);
            }
        }
    }
    printf("\n");

    finalize(cumu, pd);
}

