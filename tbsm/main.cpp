#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <chrono>
#include <iostream>
#include <fstream>
#include "bts.h"
#include "../util.h"
#include "../vector_type.h"

using namespace std;


using double9_t = vector_t<double, 9>;
using double2_t = vector_t<double, 2>;

template<size_t N>
class BtsHandler : public BtsReader {
public:
    float *curr, *next;
    int step = - 1;
    float ly, lz;

    void print_info() {
        BtsReader::print_info();
    }

    BtsHandler(string path) : BtsReader(path) {
        if (N > nv) {
            printf("BTS HANDLER ERROR: REQUIRING TOO MANY VARS!(%d > %d)\n", N, nv);
        }
        curr = new float[ny*nz*nv];
        next = new float[ny*nz*nv];
        ly = (ny - 1)*dy;
        lz = (nz - 1)*dz;
    }

    ~BtsHandler() {
        delete[] curr;
        delete[] next;
    }

    void apply_inflow_with_interpolation(
        vector_t<float, N> *inflow,
        float dyin,
        float dzin, 
        float oyin,
        float ozin,
        int nyin,
        int nzin,
        double t,
        double u_scale = 1.
    ) {
        if (step != int(t/dt)) {
            step = int(t/dt);
            load_to_array(curr, step);
            load_to_array(next, step + 1);
        }
        for (int j = 0; j < nyin; j ++) {
        for (int k = 0; k < nzin; k ++) {
            float y = j*dyin + oyin;
            float z = k*dzin + ozin;
            int jfloor = int(y/dy);
            int kfloor = int(z/dz);
            for (int v = 0; v < N; v ++) {
                float v0 = curr[btsid(jfloor, kfloor, v)];
                float v1 = curr[btsid(jfloor + 1, kfloor, v)];
                float v2 = curr[btsid(jfloor, kfloor + 1, v)];
                float v3 = curr[btsid(jfloor + 1, kfloor + 1, v)];
                float v_curr = bilinear_interpolate(
                    v0, v1, v2, v3,
                    jfloor*dy, (jfloor + 1)*dy, kfloor*dz, (kfloor + 1)*dz,
                    y, z
                );
                v0 = next[btsid(jfloor, kfloor, v)];
                v1 = next[btsid(jfloor + 1, kfloor, v)];
                v2 = next[btsid(jfloor, kfloor + 1, v)];
                v3 = next[btsid(jfloor + 1, kfloor + 1, v)];
                float v_next = bilinear_interpolate(
                    v0, v1, v2, v3,
                    jfloor*dy, (jfloor + 1)*dy, kfloor*dz, (kfloor + 1)*dz,
                    y, z
                );
                float at = (t - step*dt)/dt;
                float v_now = (at*v_next + (1 - at)*v_curr)/u_scale;
                inflow[j*nzin + k][v] = v_now;
            }
        }}
    }
};

template<size_t N>
class InflowHandler : public BtsHandler<N> {
public:
    vector_t<float, N> *inflow;
    int ny, nz;
    float dy, dz;
    float oy, oz;
    double u_scale;

    InflowHandler(int ny, int nz, string path) : BtsHandler<N>(path), ny(ny), nz(nz), dy(dy), dz(dz){
        dy = BtsHandler<N>::ly/ny;
        dz = BtsHandler<N>::lz/nz;
        float ly = (ny - 1)*dy;
        float lz = (nz - 1)*dz;
        if (ly >= BtsHandler<N>::ly || lz >= BtsHandler<N>::lz) {
            printf("INFLOW SPACE RANGE ERROR: (%f %f) >= (%f %f)\n", ly, lz, BtsHandler<N>::ly, BtsHandler<N>::lz);
        }

        oy = 0.5*(BtsHandler<N>::ly - ly);
        oz = 0.5*(BtsHandler<N>::lz - lz);

        inflow = new vector_t<float, N>[ny*nz];

        #pragma acc enter data \
        copyin(this[0:1]) create(inflow[:ny*nz])
    }

    ~InflowHandler() {
        delete[] inflow;

        #pragma acc exit data \
        delete(inflow[:ny*nz], this[0:1])
    }

    void apply_inflow(float t) {
        BtsHandler<N>::apply_inflow_with_interpolation(
            inflow,
            dy, dz,
            oy, oz,
            ny, nz,
            t,
            u_scale
        );
        #pragma acc update \
        device(inflow[:ny*nz])
    }
};

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
#pragma acc declare \
copyin(Ve)

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
#pragma acc declare \
copyin(Wgt)

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
const double Cl_ = dx_/dx;
const double Ox = Ox_/Cl_;
const double Oy = Oy_/Cl_;
const double U = 0.02;
const double Cu_ = U_/U;
const double Ct_ = Cl_/Cu_;
const double dt = 1;
const double dt_ = dt*Ct_;
const double Cnu_ = Cl_*Cu_;
const double nu = nu_/Cnu_;
const double csq = 1./3.;
const double csqi = 3;
const double tau = nu*csqi + 0.5;
const double omega = 1/tau;
const double Re_lattice = U*dx/nu;
const double CFL = U_*dt_/dx_;
const double Cpressure_ = 1.*Cu_*Cu_;

const double CRC_ = 13;
const double Turbine_x_[] = {0., 5.};
const double Turbine_y_[] = {0., 1.};
const double Diameter_ = 1.5;
const double Thick_ = 0.03;

const double Ccrc_ = 1./Cl_;
const double CRC = CRC_/Ccrc_;
const double Diameter = Diameter_/Cl_;
const double Thick = Thick_/Cl_;

const double T_ = 300;
const int Nt = T_/dt_;

template<typename T>
T square(T a) {
    return a*a;
}

double intersection_length(double x0, double x1, double x2, double x3) {
    if (x1 > x2 && x0 < x3) {
        return fmin(x3-x2, fmin(x1-x0, fmin(x1-x2, x3-x0)));
    } else {
        return 0;
    }
}

class Cumu {
public:
    double9_t *f, *fpost, *fprev, *c, *cpost;
    double2_t *shift;
    int nx, ny;
    double omega;

    Cumu(int nx, int ny, double omega) : nx(nx), ny(ny), omega(omega) {
        f = new double9_t[nx*ny];
        fpost = new double9_t[nx*ny];
        fprev = new double9_t[nx*ny];
        c = new double9_t[nx*ny];
        cpost = new double9_t[nx*ny];
        shift = new double2_t[nx*ny];

        #pragma acc enter data \
        copyin(this[0:1], f[:nx*ny], fpost[:nx*ny], fprev[:nx*ny], c[:nx*ny], cpost[:nx*ny], shift[:nx*ny])
    }

    ~Cumu() {
        delete[] f;
        delete[] fpost;
        delete[] fprev;
        delete[] c;
        delete[] cpost;
        delete[] shift;

        #pragma acc exit data \
        delete(f[:nx*ny], fpost[:nx*ny], fprev[:nx*ny], c[:nx*ny], cpost[:nx*ny], shift[:nx*ny], this[0:1])
    }

    void print_info() {
        printf("CUMULANT LBM\n");
        printf("\tdomain size = (%d %d)\n", nx, ny);
        printf("\tghost cell = %d\n", Ghost_cell);
        printf("\trelaxation rate = %lf\n", omega);
        printf("\tRe = %lf\n", Re);
    }
};

class PorousDisk {
public:
    double *dfunc;
    int nx, ny;

    PorousDisk(int nx, int ny) : nx(nx), ny(ny) {
        dfunc = new double[nx*ny]{};
        for (int i = 0; i < nx; i ++) {
        for (int j = 0; j < ny; j ++) {
            double x0 = Ox + dx*(i - 1);
            double x1 = Ox + dx*(i);
            double y0 = Oy + dx*(j - 1);
            double y1 = Oy + dx*(j);
            for (int t = 0; t < sizeof(Turbine_x_)/sizeof(double); t ++) {
                double x2 = Turbine_x_[t]/Cl_ - 0.5*Thick;
                double x3 = Turbine_x_[t]/Cl_ + 0.5*Thick;
                double y2 = Turbine_y_[t]/Cl_ - 0.5*Diameter;
                double y3 = Turbine_y_[t]/Cl_ + 0.5*Diameter;
                double cover = intersection_length(x0, x1, x2, x3)*intersection_length(y0, y1, y2, y3);
                if (cover > 0) {
                    double sigma = Diameter/6;
                    double d = fabs(0.5*(y0 + y1) - Turbine_y_[t]/Cl_);
                    const int l = i*ny + j;
                    dfunc[l] = CRC*exp(- 0.5*square(d/sigma))*cover/(dx*dx);
                    break;
                }
            }
        }}

        #pragma acc enter data \
        copyin(this[0:1], dfunc[:nx*ny])
    }

    ~PorousDisk() {
        delete[] dfunc;

        #pragma acc exit data \
        delete(dfunc[:nx*ny], this[0:1])
    }

    void print_info() {
        printf("PD INFO\n");
        printf("\tGlobal size = (%d %d)\n", nx, ny);
        printf("\tCRC = %lf\n", CRC_);
    }

    void output_dfunc(std::string path) {
        FILE *file = fopen(path.c_str(), "w");
        fprintf(file, "x,y,z,dfunc\n");
        for (int j = 1; j < ny - 1; j ++) {
        for (int i = 1; i < nx - 1; i ++) {
            fprintf(
                file,
                "%e,%e,%e,%e\n",
                (i - 0.5)*dx_ + Ox_, (j - 0.5)*dx_ + Oy_, 0.0,
                dfunc[i*ny + j]*Ccrc_
            );
        }}
        fclose(file);
    }
};

void get_macroscopic_val(
    const double9_t f,
    const double2_t shift,
    double &u,
    double &v,
    double &density
) {
    double m10 = 0, m01 = 0, m00 = 0;
    for (int q = 0; q < 9; q ++) {
        m10 += Ve[q][0]*f[q];
        m01 += Ve[q][1]*f[q];
        m00 += f[q];
    }
    u = (m10 +shift[0])/m00;
    v = (m01 +shift[1])/m00;
    density = m00;
}

double9_t get_eq_cumulant(
    const double u,
    const double v,
    const double density
) {
    return double9_t{{
        /** 00       01          02         10      11 12      20      21 22 */
        density, density*v, csq*density, density*u, 0, 0, csq*density, 0, 0
    }};
}

double9_t pdf_to_raw(const double9_t f) {
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

double9_t raw_to_pdf(const double9_t m) {
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

double9_t raw_to_central(const double9_t m, const double u, const double v) {
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

double9_t central_to_raw(const double9_t k, const double u, const double v) {
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

double9_t central_to_cumulant(const double9_t k, const double u, const double v) {
    double9_t c;
    const double density = k[0];
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

double9_t cumulant_to_central(const double9_t c, const double2_t shift) {
    double9_t k;
    const double density = c[0];
    k[0] = density;
    k[1] = shift[1];
    k[3] = shift[0];
    k[2] = c[2];
    k[4] = c[4];
    k[6] = c[6];
    k[5] = c[5];
    k[7] = c[7];
    k[8] = c[8] + (2*k[4]*k[4] + k[6]*k[2])/density;
    return k;
}

double9_t pdf_to_cumulant(
    const double9_t f,
    const double2_t shift
) {
    double u, v, density;
    get_macroscopic_val(f, shift, u, v, density);
    return central_to_cumulant(raw_to_central(pdf_to_raw(f), u, v), u, v);
}

double9_t cumulant_to_pdf(
    const double9_t c,
    const double2_t shift
) {
    double u = c[3]/c[0], v = c[1]/c[0];
    return raw_to_pdf(central_to_raw(cumulant_to_central(c, shift), u, v));
}

void compute_cumulant(
    const double9_t *f,
    const double2_t *shift,
    double9_t *c,
    const int nx,
    const int ny
) {
    #pragma acc parallel loop independent \
    present(f[:nx*ny], c[:nx*ny], shift[:nx*ny], Ve) \
    firstprivate(nx, ny)
    for (int i = 0; i < nx*ny; i ++) {
        c[i] = pdf_to_cumulant(f[i], shift[i]);
    }
}

double9_t relax_cumulant(
    const double9_t c,
    const double omega
) {
    double9_t cp;
    const double density = c[0];
    cp[0] = c[0];
    cp[1] = c[1];
    cp[3] = c[3];
    cp[2] = 0.5*(1 - omega)*(c[2] - c[6]) + density*csq;
    cp[4] = (1 - omega)*c[4];
    cp[6] = 0.5*(1 - omega)*(c[6] - c[2]) + density*csq;
    cp[5] = 0;
    cp[7] = 0;
    cp[8] = 0;
    return cp;
}

void do_collision(
    const double9_t *c,
    const double2_t *shift,
    double9_t *cp,
    const double omega,
    const int nx,
    const int ny
) {
    #pragma acc parallel loop independent \
    present(c[:nx*ny], shift[:nx*ny], cp[:nx*ny]) \
    firstprivate(omega, nx, ny)
    for (int i = 0; i < nx*ny; i ++) {
        cp[i] = relax_cumulant(c[i], omega);
    }
}

void compute_post_pdfs(
    const double9_t *cp,
    const double2_t *shift,
    double9_t *fp,
    const int nx,
    const int ny
) {
    #pragma acc parallel loop independent \
    present(cp[:nx*ny], shift[:nx*ny], fp[:nx*ny]) \
    firstprivate(nx, ny)
    for (int i = 0; i < nx*ny; i ++) {
        fp[i] = cumulant_to_pdf(cp[i], shift[i]);
    }
}

void do_streaming(
    const double9_t *fp,
    double9_t *f,
    const int nx,
    const int ny
) {
    #pragma acc parallel loop independent collapse(3) \
    present(fp[:nx*ny], f[:nx*ny]) \
    firstprivate(nx, ny)
    for (int i = 1; i < nx - 1; i ++) {
    for (int j = 1; j < ny - 1; j ++) {
    for (int q = 0; q < 9        ; q ++) {
        int src = (i - Ve[q][0])*ny + (j - Ve[q][1]);
        int dst = i*ny + j;
        f[dst][q] = fp[src][q];
    }}}
}

template<typename T>
void cpy_array(T *dst, T *src, int cnt) {
    #pragma acc parallel loop \
    present(dst[:cnt], src[:cnt]) \
    firstprivate(cnt)
    for (int i = 0; i < cnt; i ++) {
        dst[i] = src[i];
    }
}

void apply_bc(
    double9_t *f,
    const double9_t *fpost,
    const double9_t *fprev,
    const double2_t *shift,
    InflowHandler<2> *inflow_handler,
    const double t_realworld,
    const int nx,
    const int ny
) {
    /** bottom slip */
    #pragma acc parallel loop independent \
    present(f[:nx*ny], fpost[:nx*ny], Ve) \
    firstprivate(nx, ny)
    for (int i = 1; i < nx - 1; i ++) {
        const int j = 1;
        const int qlist[]{2, 5, 8};
        for (int q : qlist) {
            const int p = q - 2;
            const int dst = i*ny + j;
            const int src = (i - Ve[p][0])*ny + j;
            f[dst][q] = fpost[src][p];
        }
    }

    /** top slip */
    #pragma acc parallel loop independent \
    present(f[:nx*ny], fpost[:nx*ny], Ve) \
    firstprivate(nx, ny)
    for (int i = 1; i < nx - 1; i ++) {
        const int j = ny - 2;
        const int qlist[]{0, 3, 6};
        for (int q : qlist) {
            const int p = q + 2;
            const int dst = i*ny + j;
            const int src = (i - Ve[p][0])*ny + j;
            f[dst][q] = fpost[src][p];
        }
    }

    /** right convective outlet */
    #pragma acc parallel loop independent \
    present(f[:nx*ny], fprev[:nx*ny], shift[:nx*ny], Ve, Wgt) \
    firstprivate(nx, ny)
    for (int j = 1; j < ny - 1; j ++) {
        const int i = nx - 2;
        const int l0 = (i    )*ny + j;
        const int l1 = (i - 1)*ny + j;
        const int l2 = (i - 2)*ny + j;
        double u0, u1, u2, v0, v1, v2, density0, density1, density2;
        get_macroscopic_val(fprev[l0], shift[l0], u0, v0, density0);
        get_macroscopic_val(fprev[l1], shift[l1], u1, v1, density1);
        get_macroscopic_val(fprev[l2], shift[l2], u2, v2, density2);
        double2_t gradient{{
            3*u0*density0 - 4*u1*density1 + u2*density2,
            3*v0*density0 - 4*v1*density1 + v2*density2
        }};
        const int qlist[]{0, 1, 2};
        for (int q : qlist) {
            f[l0][q] = fprev[l0][q] - 1.5*Wgt[q]*u0*(gradient[0]*Ve[q][0] + gradient[1]*Ve[q][1]);
        }
    }

    /** left time varying velocity inlet */
    inflow_handler->apply_inflow(t_realworld);
    vector_t<float, 2> *inflow = inflow_handler->inflow;
    #pragma acc parallel loop independent \
    present(f[:nx*ny], fpost[:nx*ny], Wgt, Ve, inflow[:(ny - 2)]) \
    firstprivate(nx, ny)
    for (int j = 1; j < ny - 1; j ++) {
        const int i = 1;
        const int l = i*ny + j;
        double u, v, density;
        get_macroscopic_val(f[l], shift[l], u, v, density);
        double inflow_u = inflow[j + 1][0];
        double inflow_v = inflow[j + 1][1];
        const int qlist[]{6, 7, 8};
        for (int q : qlist) {
            const int lnk = 8 - q;
            f[l][q] = fpost[l][lnk] - 2*Wgt[lnk]*density*(Ve[lnk][0]*inflow_u + Ve[lnk][1]*inflow_v)*csqi;
        }
    }
}

void compute_shift(
    const double9_t *f,
    const double *dfunc,
    double2_t *shift,
    const int nx,
    const int ny
) {
    #pragma acc parallel loop independent \
    present(f[:nx*ny], dfunc[:nx*ny], shift[:nx*ny]) \
    firstprivate(nx, ny)
    for (int i = 0; i < nx*ny; i ++) {
        double u, v, density;
        get_macroscopic_val(f[i], shift[i], u, v, density);
        double Uabs = sqrt(u*u + v*v);
        shift[i][0] = - 0.5*u*Uabs*dfunc[i];
        shift[i][1] = - 0.5*v*Uabs*dfunc[i];
    }
}

void output(Cumu *cumu, std::string path) {
    int nx = cumu->nx, ny = cumu->ny;
    #pragma acc update \
    self(cumu->f[:nx*ny], cumu->shift[:nx*ny])
    
    FILE *file = fopen(path.c_str(), "w");
    fprintf(file, "x,y,z,u,v,w,p\n");
    for (int j = 1; j < ny - 1; j ++) {
    for (int i = 1; i < nx - 1; i ++) {
        const int l = i *ny + j;
        double u, v, density;
        get_macroscopic_val(cumu->f[l], cumu->shift[l], u, v, density);
        double p = (density - 1.)*csq;
        fprintf(
            file,
            "%e,%e,%e,%lf,%lf,%lf,%lf\n",
            (i - 0.5)*dx_ + Ox_, (j - 0.5)*dx_ + Oy_, 0.,
            u*Cu_, v*Cu_, 0., p*Cpressure_
        );
    }}
    fclose(file);
}

void init(int nx, int ny, double tau, Cumu **pcumu, PorousDisk **ppd, InflowHandler<2> **pinflow) {
    InflowHandler<2> *inflow = new InflowHandler<2>(ny - 2, 1, "TurbSim.bts");
    inflow->print_info();

    PorousDisk *pd = new PorousDisk(nx, ny);
    pd->print_info();
    pd->output_dfunc("data/dfunc.csv");

    Cumu *cumu = new Cumu(nx, ny, 1./tau);
    cumu->print_info();

    #pragma acc parallel loop independent \
    present(cumu[0:1], cumu->shift[:nx*ny], cumu->f[:nx*ny]) \
    firstprivate(nx, ny)
    for (int i = 0; i < nx*ny; i ++) {
        double2_t shift{{0., 0.}};
        double density = 1.;
        double9_t ceq = get_eq_cumulant(0., 0., density);
        cumu->shift[i] = shift;
        cumu->f[i] = cumulant_to_pdf(ceq, shift);
    }

    *pcumu = cumu;
    *ppd = pd;
    *pinflow = inflow;
}

void finalize(Cumu *cumu, PorousDisk *pd, InflowHandler<2> *inflow) {
    delete cumu;
    delete pd;
    delete inflow;
}

void main_loop(Cumu *cumu, PorousDisk *pd, InflowHandler<2> *inflow) {
    
}