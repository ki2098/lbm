#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

template <typename T>
T sq(T a) {
    return a * a;
}

template <typename T>
void copy_array(const T src[], T dst[], const int len) {
    for (int i = 0; i < len; i++) {
        dst[i] = src[i];
    }
}

template <typename T, int N>
void copy_array(const T src[][N], T dst[][N], const int len) {
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < N; j++) {
            dst[i][j] = src[i][j];
        }
    }
}

const int _LLL = 0;
const int _LLO = 1;
const int _LLR = 2;
const int _LOL = 3;
const int _LOO = 4;
const int _LOR = 5;
const int _LRL = 6;
const int _LRO = 7;
const int _LRR = 8;
const int _OLL = 9;
const int _OLO = 10;
const int _OLR = 11;
const int _OOL = 12;
const int _OOO = 13;
const int _OOR = 14;
const int _ORL = 15;
const int _ORO = 16;
const int _ORR = 17;
const int _RLL = 18;
const int _RLO = 19;
const int _RLR = 20;
const int _ROL = 21;
const int _ROO = 22;
const int _ROR = 23;
const int _RRL = 24;
const int _RRO = 25;
const int _RRR = 26;

const int _000 = 0;
const int _001 = 1;
const int _002 = 2;
const int _010 = 3;
const int _011 = 4;
const int _012 = 5;
const int _020 = 6;
const int _021 = 7;
const int _022 = 8;
const int _100 = 9;
const int _101 = 10;
const int _102 = 11;
const int _110 = 12;
const int _111 = 13;
const int _112 = 14;
const int _120 = 15;
const int _121 = 16;
const int _122 = 17;
const int _200 = 18;
const int _201 = 19;
const int _202 = 20;
const int _210 = 21;
const int _211 = 22;
const int _212 = 23;
const int _220 = 24;
const int _221 = 25;
const int _222 = 26;

const int D = 3;
const int Q = 27;

const int gc = 1;

const double L_ = 1;
const double U_ = 1;
const double U = 0.1;

const double csqi = 3;
const double csq = 1. / csqi;

struct Cumu {
    double (*f)[Q], (*ft)[Q], (*c)[Q], (*ct)[Q], (*shift)[D];
    int size[D];
    double omega;

    Cumu(const int size[D], const double omega) {
        int len = size[0]*size[1]*size[2];
        copy_array(size, this->size, D);
        this->omega = omega;
        f = new double[len][Q]();
        ft = new double[len][Q]();
        c = new double[len][Q]();
        ct = new double[len][Q]();
    }

    ~Cumu() {
        delete[] f;
        delete[] ft;
        delete[] c;
        delete[] ct;
    }

    void print_info() {
        printf("CUMULANT LBM\n");
        printf("\tdomain size = (%d %d %d)\n", size[0], size[1], size[2]);
        printf("\tguide cell = %d\n", gc);
        printf("\trelaxation rate = %ld\n", omega);
    }
};

const int Vel[Q][D] = {
    {-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0}, {-1, 0, 1}, {-1, 1, -1}, {-1, 1, 0}, {-1, 1, 1},
    {0, -1, -1},  {0, -1, 0},  {0, -1, 1},  {0, 0, -1},  {0, 0, 0},  {0, 0, 1},  {0, 1, -1},  {0, 1, 0},  {0, 1, 1},
    {1, -1, -1},  {1, -1, 0},  {1, -1, 1},  {1, 0, -1},  {1, 0, 0},  {1, 0, 1},  {1, 1, -1},  {1, 1, 0},  {1, 1, 1},
};

const double Wght[Q] = {
    1. / 216., 1. / 54., 1. / 216., 1. / 54., 2. / 27., 1. / 54., 1. / 216., 1. / 54., 1. / 216.,
    1. / 54.,  2. / 27., 1. / 54.,  2. / 27., 8. / 27., 2. / 27., 1. / 54.,  2. / 27., 1. / 54.,
    1. / 216., 1. / 54., 1. / 216., 1. / 54., 2. / 27., 1. / 54., 1. / 216., 1. / 54., 1. / 216.,
};

int index(const int i, const int j, const int k, const int size[D]) {
    return i * size[1] * size[2] + j * size[2] + k;
}

int link(const int q) {
    return (Q - 1) - q;
}

void pdf_to_macroscopic(const double f[Q], const double shift[D], double U[D], double &rho) {
    double m000 = 0;
    double m100 = 0;
    double m010 = 0;
    double m001 = 0;
    for (int q = 0; q < Q; q++) {
        m000 += f[q];
        m100 += f[q] * Vel[q][0];
        m010 += f[q] * Vel[q][1];
        m001 += f[q] * Vel[q][2];
    }
    rho = m000;
    U[0] = (m100 + shift[0]) / m000;
    U[1] = (m010 + shift[1]) / m000;
    U[2] = (m001 + shift[2]) / m000;
}

void get_eq_cumulant(const double U[D], const double rho, double c[Q]) {
    c[_000] = rho;
    c[_100] = rho * U[0];
    c[_010] = rho * U[1];
    c[_001] = rho * U[2];
    c[_200] = rho * csq;
    c[_020] = rho * csq;
    c[_002] = rho * csq;
    c[_110] = 0;
    c[_101] = 0;
    c[_011] = 0;
    c[_210] = 0;
    c[_201] = 0;
    c[_021] = 0;
    c[_120] = 0;
    c[_102] = 0;
    c[_012] = 0;
    c[_111] = 0;
    c[_211] = 0;
    c[_121] = 0;
    c[_112] = 0;
    c[_220] = 0;
    c[_202] = 0;
    c[_022] = 0;
    c[_221] = 0;
    c[_212] = 0;
    c[_122] = 0;
    c[_222] = 0;
}

void pdf_to_raw_moment(const double f[Q], double m[Q]) {
    m[_000] = f[_LLL] + f[_LLO] + f[_LLR] + f[_LOL] + f[_LOO] + f[_LOR] + f[_LRL] + f[_LRO] + f[_LRR] + f[_OLL] +
              f[_OLO] + f[_OLR] + f[_OOL] + f[_OOO] + f[_OOR] + f[_ORL] + f[_ORO] + f[_ORR] + f[_RLL] + f[_RLO] +
              f[_RLR] + f[_ROL] + f[_ROO] + f[_ROR] + f[_RRL] + f[_RRO] + f[_RRR];
    m[_001] = -f[_LLL] + f[_LLR] - f[_LOL] + f[_LOR] - f[_LRL] + f[_LRR] - f[_OLL] + f[_OLR] - f[_OOL] + f[_OOR] -
              f[_ORL] + f[_ORR] - f[_RLL] + f[_RLR] - f[_ROL] + f[_ROR] - f[_RRL] + f[_RRR];
    m[_002] = f[_LLL] + f[_LLR] + f[_LOL] + f[_LOR] + f[_LRL] + f[_LRR] + f[_OLL] + f[_OLR] + f[_OOL] + f[_OOR] +
              f[_ORL] + f[_ORR] + f[_RLL] + f[_RLR] + f[_ROL] + f[_ROR] + f[_RRL] + f[_RRR];
    m[_010] = -f[_LLL] - f[_LLO] - f[_LLR] + f[_LRL] + f[_LRO] + f[_LRR] - f[_OLL] - f[_OLO] - f[_OLR] + f[_ORL] +
              f[_ORO] + f[_ORR] - f[_RLL] - f[_RLO] - f[_RLR] + f[_RRL] + f[_RRO] + f[_RRR];
    m[_011] = f[_LLL] - f[_LLR] - f[_LRL] + f[_LRR] + f[_OLL] - f[_OLR] - f[_ORL] + f[_ORR] + f[_RLL] - f[_RLR] -
              f[_RRL] + f[_RRR];
    m[_012] = -f[_LLL] - f[_LLR] + f[_LRL] + f[_LRR] - f[_OLL] - f[_OLR] + f[_ORL] + f[_ORR] - f[_RLL] - f[_RLR] +
              f[_RRL] + f[_RRR];
    m[_020] = f[_LLL] + f[_LLO] + f[_LLR] + f[_LRL] + f[_LRO] + f[_LRR] + f[_OLL] + f[_OLO] + f[_OLR] + f[_ORL] +
              f[_ORO] + f[_ORR] + f[_RLL] + f[_RLO] + f[_RLR] + f[_RRL] + f[_RRO] + f[_RRR];
    m[_021] = -f[_LLL] + f[_LLR] - f[_LRL] + f[_LRR] - f[_OLL] + f[_OLR] - f[_ORL] + f[_ORR] - f[_RLL] + f[_RLR] -
              f[_RRL] + f[_RRR];
    m[_022] = f[_LLL] + f[_LLR] + f[_LRL] + f[_LRR] + f[_OLL] + f[_OLR] + f[_ORL] + f[_ORR] + f[_RLL] + f[_RLR] +
              f[_RRL] + f[_RRR];
    m[_100] = -f[_LLL] - f[_LLO] - f[_LLR] - f[_LOL] - f[_LOO] - f[_LOR] - f[_LRL] - f[_LRO] - f[_LRR] + f[_RLL] +
              f[_RLO] + f[_RLR] + f[_ROL] + f[_ROO] + f[_ROR] + f[_RRL] + f[_RRO] + f[_RRR];
    m[_101] = f[_LLL] - f[_LLR] + f[_LOL] - f[_LOR] + f[_LRL] - f[_LRR] - f[_RLL] + f[_RLR] - f[_ROL] + f[_ROR] -
              f[_RRL] + f[_RRR];
    m[_102] = -f[_LLL] - f[_LLR] - f[_LOL] - f[_LOR] - f[_LRL] - f[_LRR] + f[_RLL] + f[_RLR] + f[_ROL] + f[_ROR] +
              f[_RRL] + f[_RRR];
    m[_110] = f[_LLL] + f[_LLO] + f[_LLR] - f[_LRL] - f[_LRO] - f[_LRR] - f[_RLL] - f[_RLO] - f[_RLR] + f[_RRL] +
              f[_RRO] + f[_RRR];
    m[_111] = -f[_LLL] + f[_LLR] + f[_LRL] - f[_LRR] + f[_RLL] - f[_RLR] - f[_RRL] + f[_RRR];
    m[_112] = f[_LLL] + f[_LLR] - f[_LRL] - f[_LRR] - f[_RLL] - f[_RLR] + f[_RRL] + f[_RRR];
    m[_120] = -f[_LLL] - f[_LLO] - f[_LLR] - f[_LRL] - f[_LRO] - f[_LRR] + f[_RLL] + f[_RLO] + f[_RLR] + f[_RRL] +
              f[_RRO] + f[_RRR];
    m[_121] = f[_LLL] - f[_LLR] + f[_LRL] - f[_LRR] - f[_RLL] + f[_RLR] - f[_RRL] + f[_RRR];
    m[_122] = -f[_LLL] - f[_LLR] - f[_LRL] - f[_LRR] + f[_RLL] + f[_RLR] + f[_RRL] + f[_RRR];
    m[_200] = f[_LLL] + f[_LLO] + f[_LLR] + f[_LOL] + f[_LOO] + f[_LOR] + f[_LRL] + f[_LRO] + f[_LRR] + f[_RLL] +
              f[_RLO] + f[_RLR] + f[_ROL] + f[_ROO] + f[_ROR] + f[_RRL] + f[_RRO] + f[_RRR];
    m[_201] = -f[_LLL] + f[_LLR] - f[_LOL] + f[_LOR] - f[_LRL] + f[_LRR] - f[_RLL] + f[_RLR] - f[_ROL] + f[_ROR] -
              f[_RRL] + f[_RRR];
    m[_202] = f[_LLL] + f[_LLR] + f[_LOL] + f[_LOR] + f[_LRL] + f[_LRR] + f[_RLL] + f[_RLR] + f[_ROL] + f[_ROR] +
              f[_RRL] + f[_RRR];
    m[_210] = -f[_LLL] - f[_LLO] - f[_LLR] + f[_LRL] + f[_LRO] + f[_LRR] - f[_RLL] - f[_RLO] - f[_RLR] + f[_RRL] +
              f[_RRO] + f[_RRR];
    m[_211] = f[_LLL] - f[_LLR] - f[_LRL] + f[_LRR] + f[_RLL] - f[_RLR] - f[_RRL] + f[_RRR];
    m[_212] = -f[_LLL] - f[_LLR] + f[_LRL] + f[_LRR] - f[_RLL] - f[_RLR] + f[_RRL] + f[_RRR];
    m[_220] = f[_LLL] + f[_LLO] + f[_LLR] + f[_LRL] + f[_LRO] + f[_LRR] + f[_RLL] + f[_RLO] + f[_RLR] + f[_RRL] +
              f[_RRO] + f[_RRR];
    m[_221] = -f[_LLL] + f[_LLR] - f[_LRL] + f[_LRR] - f[_RLL] + f[_RLR] - f[_RRL] + f[_RRR];
    m[_222] = f[_LLL] + f[_LLR] + f[_LRL] + f[_LRR] + f[_RLL] + f[_RLR] + f[_RRL] + f[_RRR];
}

void raw_moment_to_pdf(const double m[Q], double f[Q]) {
    f[_LLL] = -0.125 * (m[_111] - m[_112] - m[_121] + m[_122] - m[_211] + m[_212] + m[_221] - m[_222]);
    f[_LLO] = 0.25 * (m[_110] - m[_112] - m[_120] + m[_122] - m[_210] + m[_212] + m[_220] - m[_222]);
    f[_LLR] = 0.125 * (m[_111] + m[_112] - m[_121] - m[_122] - m[_211] - m[_212] + m[_221] + m[_222]);
    f[_LOL] = 0.25 * (m[_101] - m[_102] - m[_121] + m[_122] - m[_201] + m[_202] + m[_221] - m[_222]);
    f[_LOO] = -0.5 * (m[_100] - m[_102] - m[_120] + m[_122] - m[_200] + m[_202] + m[_220] - m[_222]);
    f[_LOR] = -0.25 * (m[_101] + m[_102] - m[_121] - m[_122] - m[_201] - m[_202] + m[_221] + m[_222]);
    f[_LRL] = 0.125 * (m[_111] - m[_112] + m[_121] - m[_122] - m[_211] + m[_212] - m[_221] + m[_222]);
    f[_LRO] = -0.25 * (m[_110] - m[_112] + m[_120] - m[_122] - m[_210] + m[_212] - m[_220] + m[_222]);
    f[_LRR] = -0.125 * (m[_111] + m[_112] + m[_121] + m[_122] - m[_211] - m[_212] - m[_221] - m[_222]);
    f[_OLL] = 0.25 * (m[_011] - m[_012] - m[_021] + m[_022] - m[_211] + m[_212] + m[_221] - m[_222]);
    f[_OLO] = -0.5 * (m[_010] - m[_012] - m[_020] + m[_022] - m[_210] + m[_212] + m[_220] - m[_222]);
    f[_OLR] = -0.25 * (m[_011] + m[_012] - m[_021] - m[_022] - m[_211] - m[_212] + m[_221] + m[_222]);
    f[_OOL] = -0.5 * (m[_001] - m[_002] - m[_021] + m[_022] - m[_201] + m[_202] + m[_221] - m[_222]);
    f[_OOO] = 1.0 * (m[_000] - m[_002] - m[_020] + m[_022] - m[_200] + m[_202] + m[_220] - m[_222]);
    f[_OOR] = 0.5 * (m[_001] + m[_002] - m[_021] - m[_022] - m[_201] - m[_202] + m[_221] + m[_222]);
    f[_ORL] = -0.25 * (m[_011] - m[_012] + m[_021] - m[_022] - m[_211] + m[_212] - m[_221] + m[_222]);
    f[_ORO] = 0.5 * (m[_010] - m[_012] + m[_020] - m[_022] - m[_210] + m[_212] - m[_220] + m[_222]);
    f[_ORR] = 0.25 * (m[_011] + m[_012] + m[_021] + m[_022] - m[_211] - m[_212] - m[_221] - m[_222]);
    f[_RLL] = 0.125 * (m[_111] - m[_112] - m[_121] + m[_122] + m[_211] - m[_212] - m[_221] + m[_222]);
    f[_RLO] = -0.25 * (m[_110] - m[_112] - m[_120] + m[_122] + m[_210] - m[_212] - m[_220] + m[_222]);
    f[_RLR] = -0.125 * (m[_111] + m[_112] - m[_121] - m[_122] + m[_211] + m[_212] - m[_221] - m[_222]);
    f[_ROL] = -0.25 * (m[_101] - m[_102] - m[_121] + m[_122] + m[_201] - m[_202] - m[_221] + m[_222]);
    f[_ROO] = 0.5 * (m[_100] - m[_102] - m[_120] + m[_122] + m[_200] - m[_202] - m[_220] + m[_222]);
    f[_ROR] = 0.25 * (m[_101] + m[_102] - m[_121] - m[_122] + m[_201] + m[_202] - m[_221] - m[_222]);
    f[_RRL] = -0.125 * (m[_111] - m[_112] + m[_121] - m[_122] + m[_211] - m[_212] + m[_221] - m[_222]);
    f[_RRO] = 0.25 * (m[_110] - m[_112] + m[_120] - m[_122] + m[_210] - m[_212] + m[_220] - m[_222]);
    f[_RRR] = 0.125 * (m[_111] + m[_112] + m[_121] + m[_122] + m[_211] + m[_212] + m[_221] + m[_222]);
}

void raw_moment_to_central_moment(const double m[Q], const double U[D], double k[Q]) {
    double u = U[0];
    double v = U[1];
    double w = U[2];
    double u2 = u * u;
    double v2 = v * v;
    double w2 = w * w;
    k[_000] = m[_000];
    k[_001] = -m[_000] * w + m[_001];
    k[_002] = m[_000] * w2 - 2 * m[_001] * w + m[_002];
    k[_010] = -m[_000] * v + m[_010];
    k[_011] = -m[_010] * w + m[_011] + v * (m[_000] * w - m[_001]);
    k[_012] = m[_010] * w2 - 2 * m[_011] * w + m[_012] - v * (m[_000] * w2 - 2 * m[_001] * w + m[_002]);
    k[_020] = m[_000] * v2 - 2 * m[_010] * v + m[_020];
    k[_021] = -m[_020] * w + m[_021] + v2 * (-m[_000] * w + m[_001]) + 2 * v * (m[_010] * w - m[_011]);
    k[_022] = m[_020] * w2 - 2 * m[_021] * w + m[_022] + v2 * (m[_000] * w2 - 2 * m[_001] * w + m[_002]) -
              2 * v * (m[_010] * w2 - 2 * m[_011] * w + m[_012]);
    k[_100] = -m[_000] * u + m[_100];
    k[_101] = -m[_100] * w + m[_101] + u * (m[_000] * w - m[_001]);
    k[_102] = m[_100] * w2 - 2 * m[_101] * w + m[_102] - u * (m[_000] * w2 - 2 * m[_001] * w + m[_002]);
    k[_110] = -m[_100] * v + m[_110] + u * (m[_000] * v - m[_010]);
    k[_111] = -m[_110] * w + m[_111] - u * (-m[_010] * w + m[_011] + v * (m[_000] * w - m[_001])) +
              v * (m[_100] * w - m[_101]);
    k[_112] = m[_110] * w2 - 2 * m[_111] * w + m[_112] -
              u * (m[_010] * w2 - 2 * m[_011] * w + m[_012] - v * (m[_000] * w2 - 2 * m[_001] * w + m[_002])) -
              v * (m[_100] * w2 - 2 * m[_101] * w + m[_102]);
    k[_120] = m[_100] * v2 - 2 * m[_110] * v + m[_120] - u * (m[_000] * v2 - 2 * m[_010] * v + m[_020]);
    k[_121] = -m[_120] * w + m[_121] +
              u * (m[_020] * w - m[_021] + v2 * (m[_000] * w - m[_001]) - 2 * v * (m[_010] * w - m[_011])) +
              v2 * (-m[_100] * w + m[_101]) + 2 * v * (m[_110] * w - m[_111]);
    k[_122] = m[_120] * w2 - 2 * m[_121] * w + m[_122] -
              u * (m[_020] * w2 - 2 * m[_021] * w + m[_022] + v2 * (m[_000] * w2 - 2 * m[_001] * w + m[_002]) -
                   2 * v * (m[_010] * w2 - 2 * m[_011] * w + m[_012])) +
              v2 * (m[_100] * w2 - 2 * m[_101] * w + m[_102]) - 2 * v * (m[_110] * w2 - 2 * m[_111] * w + m[_112]);
    k[_200] = m[_000] * u2 - 2 * m[_100] * u + m[_200];
    k[_201] = -m[_200] * w + m[_201] + u2 * (-m[_000] * w + m[_001]) + 2 * u * (m[_100] * w - m[_101]);
    k[_202] = m[_200] * w2 - 2 * m[_201] * w + m[_202] + u2 * (m[_000] * w2 - 2 * m[_001] * w + m[_002]) -
              2 * u * (m[_100] * w2 - 2 * m[_101] * w + m[_102]);
    k[_210] = -m[_200] * v + m[_210] + u2 * (-m[_000] * v + m[_010]) + 2 * u * (m[_100] * v - m[_110]);
    k[_211] = -m[_210] * w + m[_211] + u2 * (-m[_010] * w + m[_011] + v * (m[_000] * w - m[_001])) -
              2 * u * (-m[_110] * w + m[_111] + v * (m[_100] * w - m[_101])) + v * (m[_200] * w - m[_201]);
    k[_212] = m[_210] * w2 - 2 * m[_211] * w + m[_212] +
              u2 * (m[_010] * w2 - 2 * m[_011] * w + m[_012] - v * (m[_000] * w2 - 2 * m[_001] * w + m[_002])) -
              2 * u * (m[_110] * w2 - 2 * m[_111] * w + m[_112] - v * (m[_100] * w2 - 2 * m[_101] * w + m[_102])) -
              v * (m[_200] * w2 - 2 * m[_201] * w + m[_202]);
    k[_220] = m[_200] * v2 - 2 * m[_210] * v + m[_220] + u2 * (m[_000] * v2 - 2 * m[_010] * v + m[_020]) -
              2 * u * (m[_100] * v2 - 2 * m[_110] * v + m[_120]);
    k[_221] = -m[_220] * w + m[_221] +
              u2 * (-m[_020] * w + m[_021] - v2 * (m[_000] * w - m[_001]) + 2 * v * (m[_010] * w - m[_011])) +
              2 * u * (m[_120] * w - m[_121] + v2 * (m[_100] * w - m[_101]) - 2 * v * (m[_110] * w - m[_111])) +
              v2 * (-m[_200] * w + m[_201]) + 2 * v * (m[_210] * w - m[_211]);
    k[_222] = m[_220] * w2 - 2 * m[_221] * w + m[_222] +
              u2 * (m[_020] * w2 - 2 * m[_021] * w + m[_022] + v2 * (m[_000] * w2 - 2 * m[_001] * w + m[_002]) -
                    2 * v * (m[_010] * w2 - 2 * m[_011] * w + m[_012])) -
              2 * u *
                  (m[_120] * w2 - 2 * m[_121] * w + m[_122] + v2 * (m[_100] * w2 - 2 * m[_101] * w + m[_102]) -
                   2 * v * (m[_110] * w2 - 2 * m[_111] * w + m[_112])) +
              v2 * (m[_200] * w2 - 2 * m[_201] * w + m[_202]) - 2 * v * (m[_210] * w2 - 2 * m[_211] * w + m[_212]);
}

void central_moment_to_raw_moment(const double k[Q], const double U[D], double m[Q]) {
    double u = U[0];
    double v = U[1];
    double w = U[2];
    double u2 = u * u;
    double v2 = v * v;
    double w2 = w * w;
    m[_000] = k[_000];
    m[_001] = k[_000] * w + k[_001];
    m[_002] = k[_000] * w2 + 2 * k[_001] * w + k[_002];
    m[_010] = k[_000] * v + k[_010];
    m[_011] = k[_010] * w + k[_011] + v * (k[_000] * w + k[_001]);
    m[_012] = k[_010] * w2 + 2 * k[_011] * w + k[_012] + v * (k[_000] * w2 + 2 * k[_001] * w + k[_002]);
    m[_020] = k[_000] * v2 + 2 * k[_010] * v + k[_020];
    m[_021] = k[_020] * w + k[_021] + v2 * (k[_000] * w + k[_001]) + 2 * v * (k[_010] * w + k[_011]);
    m[_022] = k[_020] * w2 + 2 * k[_021] * w + k[_022] + v2 * (k[_000] * w2 + 2 * k[_001] * w + k[_002]) +
              2 * v * (k[_010] * w2 + 2 * k[_011] * w + k[_012]);
    m[_100] = k[_000] * u + k[_100];
    m[_101] = k[_100] * w + k[_101] + u * (k[_000] * w + k[_001]);
    m[_102] = k[_100] * w2 + 2 * k[_101] * w + k[_102] + u * (k[_000] * w2 + 2 * k[_001] * w + k[_002]);
    m[_110] = k[_100] * v + k[_110] + u * (k[_000] * v + k[_010]);
    m[_111] =
        k[_110] * w + k[_111] + u * (k[_010] * w + k[_011] + v * (k[_000] * w + k[_001])) + v * (k[_100] * w + k[_101]);
    m[_112] = k[_110] * w2 + 2 * k[_111] * w + k[_112] +
              u * (k[_010] * w2 + 2 * k[_011] * w + k[_012] + v * (k[_000] * w2 + 2 * k[_001] * w + k[_002])) +
              v * (k[_100] * w2 + 2 * k[_101] * w + k[_102]);
    m[_120] = k[_100] * v2 + 2 * k[_110] * v + k[_120] + u * (k[_000] * v2 + 2 * k[_010] * v + k[_020]);
    m[_121] = k[_120] * w + k[_121] +
              u * (k[_020] * w + k[_021] + v2 * (k[_000] * w + k[_001]) + 2 * v * (k[_010] * w + k[_011])) +
              v2 * (k[_100] * w + k[_101]) + 2 * v * (k[_110] * w + k[_111]);
    m[_122] = k[_120] * w2 + 2 * k[_121] * w + k[_122] +
              u * (k[_020] * w2 + 2 * k[_021] * w + k[_022] + v2 * (k[_000] * w2 + 2 * k[_001] * w + k[_002]) +
                   2 * v * (k[_010] * w2 + 2 * k[_011] * w + k[_012])) +
              v2 * (k[_100] * w2 + 2 * k[_101] * w + k[_102]) + 2 * v * (k[_110] * w2 + 2 * k[_111] * w + k[_112]);
    m[_200] = k[_000] * u2 + 2 * k[_100] * u + k[_200];
    m[_201] = k[_200] * w + k[_201] + u2 * (k[_000] * w + k[_001]) + 2 * u * (k[_100] * w + k[_101]);
    m[_202] = k[_200] * w2 + 2 * k[_201] * w + k[_202] + u2 * (k[_000] * w2 + 2 * k[_001] * w + k[_002]) +
              2 * u * (k[_100] * w2 + 2 * k[_101] * w + k[_102]);
    m[_210] = k[_200] * v + k[_210] + u2 * (k[_000] * v + k[_010]) + 2 * u * (k[_100] * v + k[_110]);
    m[_211] = k[_210] * w + k[_211] + u2 * (k[_010] * w + k[_011] + v * (k[_000] * w + k[_001])) +
              2 * u * (k[_110] * w + k[_111] + v * (k[_100] * w + k[_101])) + v * (k[_200] * w + k[_201]);
    m[_212] = k[_210] * w2 + 2 * k[_211] * w + k[_212] +
              u2 * (k[_010] * w2 + 2 * k[_011] * w + k[_012] + v * (k[_000] * w2 + 2 * k[_001] * w + k[_002])) +
              2 * u * (k[_110] * w2 + 2 * k[_111] * w + k[_112] + v * (k[_100] * w2 + 2 * k[_101] * w + k[_102])) +
              v * (k[_200] * w2 + 2 * k[_201] * w + k[_202]);
    m[_220] = k[_200] * v2 + 2 * k[_210] * v + k[_220] + u2 * (k[_000] * v2 + 2 * k[_010] * v + k[_020]) +
              2 * u * (k[_100] * v2 + 2 * k[_110] * v + k[_120]);
    m[_221] = k[_220] * w + k[_221] +
              u2 * (k[_020] * w + k[_021] + v2 * (k[_000] * w + k[_001]) + 2 * v * (k[_010] * w + k[_011])) +
              2 * u * (k[_120] * w + k[_121] + v2 * (k[_100] * w + k[_101]) + 2 * v * (k[_110] * w + k[_111])) +
              v2 * (k[_200] * w + k[_201]) + 2 * v * (k[_210] * w + k[_211]);
    m[_222] = k[_220] * w2 + 2 * k[_221] * w + k[_222] +
              u2 * (k[_020] * w2 + 2 * k[_021] * w + k[_022] + v2 * (k[_000] * w2 + 2 * k[_001] * w + k[_002]) +
                    2 * v * (k[_010] * w2 + 2 * k[_011] * w + k[_012])) +
              2 * u *
                  (k[_120] * w2 + 2 * k[_121] * w + k[_122] + v2 * (k[_100] * w2 + 2 * k[_101] * w + k[_102]) +
                   2 * v * (k[_110] * w2 + 2 * k[_111] * w + k[_112])) +
              v2 * (k[_200] * w2 + 2 * k[_201] * w + k[_202]) + 2 * v * (k[_210] * w2 + 2 * k[_211] * w + k[_212]);
}

void central_moment_to_cumulant(const double k[Q], const double U[D], double c[Q]) {
    double rho = k[_000];
    c[_000] = rho;
    c[_100] = rho * U[0];
    c[_010] = rho * U[1];
    c[_001] = rho * U[2];
    c[_110] = k[_110];
    c[_101] = k[_101];
    c[_011] = k[_011];
    c[_200] = k[_200];
    c[_020] = k[_020];
    c[_002] = k[_002];
    c[_210] = k[_210];
    c[_201] = k[_201];
    c[_021] = k[_021];
    c[_120] = k[_120];
    c[_102] = k[_102];
    c[_012] = k[_012];
    c[_111] = k[_111];
    c[_211] = k[_211] - (k[_200] * k[_011] + 2 * k[_110] * k[_101]) / rho;
    c[_121] = k[_121] - (k[_020] * k[_101] + 2 * k[_110] * k[_011]) / rho;
    c[_112] = k[_112] - (k[_002] * k[_110] + 2 * k[_101] * k[_011]) / rho;
    c[_220] = k[_220] - (k[_200] * k[_020] + 2 * sq(k[_110])) / rho;
    c[_202] = k[_202] - (k[_200] * k[_002] + 2 * sq(k[_101])) / rho;
    c[_022] = k[_022] - (k[_020] * k[_002] + 2 * sq(k[_011])) / rho;
    c[_122] = k[_122] - (k[_002] * k[_120] + k[_020] * k[_102] + 4 * k[_011] * k[_111] +
                         2 * (k[_101] * k[_021] + k[_110] * k[_012])) /
                            rho;
    c[_212] = k[_212] - (k[_002] * k[_210] + k[_200] * k[_012] + 4 * k[_101] * k[_111] +
                         2 * (k[_011] * k[_201] + k[_110] * k[_102])) /
                            rho;
    c[_221] = k[_221] - (k[_020] * k[_201] + k[_200] * k[_021] + 4 * k[_110] * k[_111] +
                         2 * (k[_011] * k[_210] + k[_101] * k[_120])) /
                            rho;
    c[_222] = k[_222] -
              (4 * sq(k[_111]) + k[_200] * k[_022] + k[_020] * k[_202] + k[_002] * k[_220] +
               4 * (k[_011] * k[_211] + k[_101] * k[_121] + k[_110] * k[_112]) +
               2 * (k[_120] * k[_102] + k[_210] * k[_012] + k[_201] * k[_021])) /
                  rho +
              (16 * k[_110] * k[_101] * k[_011] +
               4 * (sq(k[_101]) * k[_020] + sq(k[_011]) * k[_200] + sq(k[_110]) * k[_002]) +
               2 * k[_200] * k[_020] * k[_002]) /
                  sq(rho);
}

void cumulant_to_central_moment(const double c[Q], const double shift[D], double k[Q]) {
    double rho = c[_000];
    k[_000] = rho;
    k[_100] = shift[0];
    k[_010] = shift[1];
    k[_001] = shift[2];
    k[_110] = c[_110];
    k[_101] = c[_101];
    k[_011] = c[_011];
    k[_200] = c[_200];
    k[_020] = c[_020];
    k[_002] = c[_002];
    k[_210] = c[_210];
    k[_201] = c[_201];
    k[_021] = c[_021];
    k[_120] = c[_120];
    k[_102] = c[_102];
    k[_012] = c[_012];
    k[_111] = c[_111];
    k[_211] = c[_211] + (k[_200] * k[_011] + 2 * k[_110] * k[_101]) / rho;
    k[_121] = c[_121] + (k[_020] * k[_101] + 2 * k[_110] * k[_011]) / rho;
    k[_112] = c[_112] + (k[_002] * k[_110] + 2 * k[_101] * k[_011]) / rho;
    k[_220] = c[_220] + (k[_200] * k[_020] + 2 * sq(k[_110])) / rho;
    k[_202] = c[_202] + (k[_200] * k[_002] + 2 * sq(k[_101])) / rho;
    k[_022] = c[_022] + (k[_020] * k[_002] + 2 * sq(k[_011])) / rho;
    k[_122] = c[_122] + (k[_002] * k[_120] + k[_020] * k[_102] + 4 * k[_011] * k[_111] +
                         2 * (k[_101] * k[_021] + k[_110] * k[_012])) /
                            rho;
    k[_212] = c[_212] + (k[_002] * k[_210] + k[_200] * k[_012] + 4 * k[_101] * k[_111] +
                         2 * (k[_011] * k[_201] + k[_110] * k[_102])) /
                            rho;
    k[_221] = c[_221] + (k[_020] * k[_201] + k[_200] * k[_021] + 4 * k[_110] * k[_111] +
                         2 * (k[_011] * k[_210] + k[_101] * k[_120])) /
                            rho;
    k[_222] = c[_222] +
              (4 * sq(k[_111]) + k[_200] * k[_022] + k[_020] * k[_202] + k[_002] * k[_220] +
               4 * (k[_011] * k[_211] + k[_101] * k[_121] + k[_110] * k[_112]) +
               2 * (k[_120] * k[_102] + k[_210] * k[_012] + k[_201] * k[_021])) /
                  rho +
              (16 * k[_110] * k[_101] * k[_011] +
               4 * (sq(k[_101]) * k[_020] + sq(k[_011]) * k[_200] + sq(k[_110]) * k[_002]) +
               2 * k[_200] * k[_020] * k[_002]) /
                  sq(rho);
}

void pdf_to_cumulant(const double f[Q], const double shift[D], double c[Q]) {
    double U[D], rho, m[Q], k[Q];
    pdf_to_macroscopic(f, shift, U, rho);
    pdf_to_raw_moment(f, m);
    raw_moment_to_central_moment(m, U, k);
    central_moment_to_cumulant(k, U, c);
}

void cumulant_to_pdf(const double c[Q], const double shift[D], double f[Q]) {
    double U[] = {c[_100] / c[_000], c[_010] / c[_000], c[_001] / c[_000]};
    double m[Q], k[Q];
    cumulant_to_central_moment(c, shift, k);
    central_moment_to_raw_moment(k, U, m);
    raw_moment_to_pdf(m, f);
}

void relax_cumulant(const double c[Q], const double omega, double ct[Q]) {
    double rho = c[_000];
    ct[_000] = rho;
    ct[_100] = c[_100];
    ct[_010] = c[_010];
    ct[_001] = c[_001];
    ct[_110] = (1 - omega) * c[_110];
    ct[_101] = (1 - omega) * c[_101];
    ct[_011] = (1 - omega) * c[_011];
    ct[_200] = (rho + (1 - omega) * (2 * c[_200] - c[_020] - c[_002])) / 3;
    ct[_020] = (rho + (1 - omega) * (2 * c[_020] - c[_200] - c[_002])) / 3;
    ct[_002] = (rho + (1 - omega) * (2 * c[_002] - c[_200] - c[_020])) / 3;
    ct[_210] = 0;
    ct[_201] = 0;
    ct[_021] = 0;
    ct[_120] = 0;
    ct[_102] = 0;
    ct[_012] = 0;
    ct[_111] = 0;
    ct[_220] = 0;
    ct[_202] = 0;
    ct[_022] = 0;
    ct[_211] = 0;
    ct[_121] = 0;
    ct[_112] = 0;
    ct[_221] = 0;
    ct[_212] = 0;
    ct[_122] = 0;
    ct[_222] = 0;
}

void compute_pre_collision_cumulant(const double f[][Q], const double shift[][D], double c[][Q], const int size[D]) {
    int len = size[0] * size[1] * size[2];
    for (int i = 0; i < len; i++) {
        pdf_to_cumulant(f[i], shift[i], c[i]);
    }
}

void apply_collision(const double c[][Q], const double omega, double ct[][Q], const int size[D]) {
    int len = size[0] * size[1] * size[2];
    for (int i = 0; i < len; i++) {
        relax_cumulant(c[i], omega, ct[i]);
    }
}

void compute_post_collision_pdf(const double c[][Q], const double shift[][D], double f[][Q], const int size[D]) {
    int len = size[0] * size[1] * size[2];
    for (int i = 0; i < len; i++) {
        cumulant_to_pdf(c[i], shift[i], f[i]);
    }
}

void apply_streaming(const double ft[][Q], double f[][Q], const int size[D]) {
    for (int i = gc; i < size[0] - gc; i++) {
        for (int j = gc; j < size[1] - gc; j++) {
            for (int k = gc; k < size[2] - gc; k++) {
                for (int q = 0; q < Q; q++) {
                    int sid = index(i - Vel[q][0], j - Vel[q][1], k - Vel[q][2], size);
                    int did = index(i, j, k, size);
                    f[did][q] = ft[sid][q];
                }
            }
        }
    }
}

void apply_boundary_condition(const double ft[][Q], const double ulid, double f[][Q], const int size[D]) {
    // bottom wall
    for (int i = gc; i < size[0] - gc; i++) {
        for (int j = gc; j < size[1] - gc; j++) {
            int id = index(i, j, gc, size);
            int qlist[] = {_LLR, _LOR, _LRR, _OLR, _OOR, _ORR, _RLR, _ROR, _RRR};
            for (auto q : qlist) {
                f[id][q] = ft[id][link(q)];
            }
        }
    }
    // top moving wall
    for (int i = gc; i < size[0] - gc; i++) {
        for (int j = gc; j < size[1] - gc; j++) {
            int id = index(i, j, size[2] - gc - 1, size);
            int qlist[] = {_LLL, _LOL, _LRL, _OLL, _OOL, _ORL, _RLL, _ROL, _RRL};
            for (auto q : qlist) {
                int lq = link(q);
                f[id][q] = ft[id][lq] - Vel[lq][0] * 2 * ulid * Wght[lq] * csqi;
            }
        }
    }
    // left wall
    for (int j = gc; j < size[1] - gc; j++) {
        for (int k = gc; k < size[2] - gc; k++) {
            int id = index(gc, j, k, size);
            int qlist[] = {_RLL, _RLO, _RLR, _ROL, _ROO, _ROR, _RRL, _RRO, _RRR};
            for (auto q : qlist) {
                f[id][q] = ft[id][link(q)];
            }
        }
    }
    // right wall
    for (int j = gc; j < size[1] - gc; j++) {
        for (int k = gc; k < size[2] - gc; k++) {
            int id = index(size[0] - gc - 1, j, k, size);
            int qlist[] = {_LLL, _LLO, _LLR, _LOL, _LOO, _LOR, _LRL, _LRO, _LRR};
            for (auto q : qlist) {
                f[id][q] = ft[id][link(q)];
            }
        }
    }
    // back wall
    for (int i = gc; i < size[0] - gc; i++) {
        for (int k = gc; k < size[2] - gc; k++) {
            int id = index(i, gc, k, size);
            int qlist[] = {_LRL, _LRO, _LRR, _ORL, _ORO, _ORR, _RRL, _RRO, _RRR};
            for (auto q : qlist) {
                f[id][q] = ft[id][link(q)];
            }
        }
    }
    // front wall
    for (int i = gc; i < size[0] - gc; i++) {
        for (int k = gc; k < size[2] - gc; k++) {
            int id = index(i, size[1] - gc - 1, k, size);
            int qlist[] = {_LLL, _LLO, _LLR, _OLL, _OLO, _OLR, _RLL, _RLO, _RLR};
            for (auto q : qlist) {
                f[id][q] = ft[id][link(q)];
            }
        }
    }
}
