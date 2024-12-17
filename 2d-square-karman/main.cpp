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
const double Re = 10000;
const double Lx_ = 15*L_;
const double Ly_ =  5*L_;
const int Cells_per_length = 50;
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

const double T_ = 100;
const int N_step = T_/dt_;

class Cumu {
public:
    double9_t *f, *fpost, *c, *cpost;
    int imax, jmax;
    double omega;

    Cumu(int imax, int jmax, double omega) : imax(imax), jmax(jmax), omega(omega) {
        f = new double9_t[imax*jmax];
        fpost = new double9_t[imax*jmax];
        c = new double9_t[imax*jmax];
        cpost = new double9_t[imax*jmax];

        #pragma acc enter data \
        copyin(this[0:1], f[:imax*jmax], fpost[:imax*jmax], c[:imax*jmax], cpost[:imax*jmax])
    }

    ~Cumu() {
        delete[] f;
        delete[] fpost;
        delete[] c;
        delete[] cpost;

        #pragma acc exit data \
        delete(f[:imax*jmax], fpost[:imax*jmax], c[:imax*jmax], cpost[:imax*jmax], this[0:1])
    }
};



