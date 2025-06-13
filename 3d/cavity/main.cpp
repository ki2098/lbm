#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace std;

template <typename T>
T sq(T a) {
    return a * a;
}

int index(int i, int j, int k, int *size) {
    return i * size[1] * size[2] + j * size[2] + k;
}

const int Ve[Q][D] = {
    {-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0}, {-1, 0, 1}, {-1, 1, -1}, {-1, 1, 0}, {-1, 1, 1},
    {0, -1, -1},  {0, -1, 0},  {0, -1, 1},  {0, 0, -1},  {0, 0, 0},  {0, 0, 1},  {0, 1, -1},  {0, 1, 0},  {0, 1, 1},
    {1, -1, -1},  {1, -1, 0},  {1, -1, 1},  {1, 0, -1},  {1, 0, 0},  {1, 0, 1},  {1, 1, -1},  {1, 1, 0},  {1, 1, 1},
};

const double Wght[Q] = {
    1. / 216., 1. / 54., 1. / 216., 1. / 54., 2. / 27., 1. / 54., 1. / 216., 1. / 54., 1. / 216.,
    1. / 54.,  2. / 27., 1. / 54.,  2. / 27., 8. / 27., 2. / 27., 1. / 54.,  2. / 27., 1. / 54.,
    1. / 216., 1. / 54., 1. / 216., 1. / 54., 2. / 27., 1. / 54., 1. / 216., 1. / 54., 1. / 216.,
};

const int Q = 27;
const int D = 3;

int get_link_q(int q) {
    return Q - 1 - q;
}

const double CSQI = 3;

const int _lll = 0;
const int _llo = 1;
const int _llr = 2;
const int _lol = 3;
const int _loo = 4;
const int _lor = 5;
const int _lrl = 6;
const int _lro = 7;
const int _lrr = 8;
const int _oll = 9;
const int _olo = 10;
const int _olr = 11;
const int _ool = 12;
const int _ooo = 13;
const int _oor = 14;
const int _orl = 15;
const int _oro = 16;
const int _orr = 17;
const int _rll = 18;
const int _rlo = 19;
const int _rlr = 20;
const int _rol = 21;
const int _roo = 22;
const int _ror = 23;
const int _rrl = 24;
const int _rro = 25;
const int _rrr = 26;

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

void pdf_to_macroscopics(double f[Q], double U[D], double &rho) {
    double m000 = 0, m100 = 0, m010 = 0, m001 = 0;
    for (int q = 0; q < Q; q++) {
        m000 += f[q];
        m100 += f[q] * Ve[q][0];
        m010 += f[q] * Ve[q][1];
        m001 += f[q] * Ve[q][2];
    }
    U[0] = m100 / m000;
    U[1] = m010 / m000;
    U[2] = m001 / m000;
    rho = m000;
}

void get_eq_cumulant(double rho, double U[D], double ceq[Q]) {
    ceq[_000] = 0;
    ceq[_100] = rho * U[0];
    ceq[_010] = rho * U[1];
    ceq[_001] = rho * U[2];
    ceq[_200] = CSQI * rho;
    ceq[_020] = CSQI * rho;
    ceq[_002] = CSQI * rho;
    ceq[_110] = 0;
    ceq[_101] = 0;
    ceq[_011] = 0;
    ceq[_120] = 0;
    ceq[_102] = 0;
    ceq[_012] = 0;
    ceq[_210] = 0;
    ceq[_201] = 0;
    ceq[_021] = 0;
    ceq[_111] = 0;
    ceq[_220] = 0;
    ceq[_202] = 0;
    ceq[_022] = 0;
    ceq[_211] = 0;
    ceq[_121] = 0;
    ceq[_112] = 0;
    ceq[_221] = 0;
    ceq[_212] = 0;
    ceq[_122] = 0;
    ceq[_222] = 0;
}

void pdf_to_raw_moment(double f[Q], double m[Q]) {
    m[_000] = f[_lll] + f[_llo] + f[_llr] + f[_lol] + f[_loo] + f[_lor] + f[_lrl] + f[_lro] + f[_lrr] + f[_oll] +
              f[_olo] + f[_olr] + f[_ool] + f[_ooo] + f[_oor] + f[_orl] + f[_oro] + f[_orr] + f[_rll] + f[_rlo] +
              f[_rlr] + f[_rol] + f[_roo] + f[_ror] + f[_rrl] + f[_rro] + f[_rrr];
    m[_001] = -f[_lll] + f[_llr] - f[_lol] + f[_lor] - f[_lrl] + f[_lrr] - f[_oll] + f[_olr] - f[_ool] + f[_oor] -
              f[_orl] + f[_orr] - f[_rll] + f[_rlr] - f[_rol] + f[_ror] - f[_rrl] + f[_rrr];
    m[_002] = f[_lll] + f[_llr] + f[_lol] + f[_lor] + f[_lrl] + f[_lrr] + f[_oll] + f[_olr] + f[_ool] + f[_oor] +
              f[_orl] + f[_orr] + f[_rll] + f[_rlr] + f[_rol] + f[_ror] + f[_rrl] + f[_rrr];
    m[_010] = -f[_lll] - f[_llo] - f[_llr] + f[_lrl] + f[_lro] + f[_lrr] - f[_oll] - f[_olo] - f[_olr] + f[_orl] +
              f[_oro] + f[_orr] - f[_rll] - f[_rlo] - f[_rlr] + f[_rrl] + f[_rro] + f[_rrr];
    m[_011] = f[_lll] - f[_llr] - f[_lrl] + f[_lrr] + f[_oll] - f[_olr] - f[_orl] + f[_orr] + f[_rll] - f[_rlr] -
              f[_rrl] + f[_rrr];
    m[_012] = -f[_lll] - f[_llr] + f[_lrl] + f[_lrr] - f[_oll] - f[_olr] + f[_orl] + f[_orr] - f[_rll] - f[_rlr] +
              f[_rrl] + f[_rrr];
    m[_020] = f[_lll] + f[_llo] + f[_llr] + f[_lrl] + f[_lro] + f[_lrr] + f[_oll] + f[_olo] + f[_olr] + f[_orl] +
              f[_oro] + f[_orr] + f[_rll] + f[_rlo] + f[_rlr] + f[_rrl] + f[_rro] + f[_rrr];
    m[_021] = -f[_lll] + f[_llr] - f[_lrl] + f[_lrr] - f[_oll] + f[_olr] - f[_orl] + f[_orr] - f[_rll] + f[_rlr] -
              f[_rrl] + f[_rrr];
    m[_022] = f[_lll] + f[_llr] + f[_lrl] + f[_lrr] + f[_oll] + f[_olr] + f[_orl] + f[_orr] + f[_rll] + f[_rlr] +
              f[_rrl] + f[_rrr];
    m[_100] = -f[_lll] - f[_llo] - f[_llr] - f[_lol] - f[_loo] - f[_lor] - f[_lrl] - f[_lro] - f[_lrr] + f[_rll] +
              f[_rlo] + f[_rlr] + f[_rol] + f[_roo] + f[_ror] + f[_rrl] + f[_rro] + f[_rrr];
    m[_101] = f[_lll] - f[_llr] + f[_lol] - f[_lor] + f[_lrl] - f[_lrr] - f[_rll] + f[_rlr] - f[_rol] + f[_ror] -
              f[_rrl] + f[_rrr];
    m[_102] = -f[_lll] - f[_llr] - f[_lol] - f[_lor] - f[_lrl] - f[_lrr] + f[_rll] + f[_rlr] + f[_rol] + f[_ror] +
              f[_rrl] + f[_rrr];
    m[_110] = f[_lll] + f[_llo] + f[_llr] - f[_lrl] - f[_lro] - f[_lrr] - f[_rll] - f[_rlo] - f[_rlr] + f[_rrl] +
              f[_rro] + f[_rrr];
    m[_111] = -f[_lll] + f[_llr] + f[_lrl] - f[_lrr] + f[_rll] - f[_rlr] - f[_rrl] + f[_rrr];
    m[_112] = f[_lll] + f[_llr] - f[_lrl] - f[_lrr] - f[_rll] - f[_rlr] + f[_rrl] + f[_rrr];
    m[_120] = -f[_lll] - f[_llo] - f[_llr] - f[_lrl] - f[_lro] - f[_lrr] + f[_rll] + f[_rlo] + f[_rlr] + f[_rrl] +
              f[_rro] + f[_rrr];
    m[_121] = f[_lll] - f[_llr] + f[_lrl] - f[_lrr] - f[_rll] + f[_rlr] - f[_rrl] + f[_rrr];
    m[_122] = -f[_lll] - f[_llr] - f[_lrl] - f[_lrr] + f[_rll] + f[_rlr] + f[_rrl] + f[_rrr];
    m[_200] = f[_lll] + f[_llo] + f[_llr] + f[_lol] + f[_loo] + f[_lor] + f[_lrl] + f[_lro] + f[_lrr] + f[_rll] +
              f[_rlo] + f[_rlr] + f[_rol] + f[_roo] + f[_ror] + f[_rrl] + f[_rro] + f[_rrr];
    m[_201] = -f[_lll] + f[_llr] - f[_lol] + f[_lor] - f[_lrl] + f[_lrr] - f[_rll] + f[_rlr] - f[_rol] + f[_ror] -
              f[_rrl] + f[_rrr];
    m[_202] = f[_lll] + f[_llr] + f[_lol] + f[_lor] + f[_lrl] + f[_lrr] + f[_rll] + f[_rlr] + f[_rol] + f[_ror] +
              f[_rrl] + f[_rrr];
    m[_210] = -f[_lll] - f[_llo] - f[_llr] + f[_lrl] + f[_lro] + f[_lrr] - f[_rll] - f[_rlo] - f[_rlr] + f[_rrl] +
              f[_rro] + f[_rrr];
    m[_211] = f[_lll] - f[_llr] - f[_lrl] + f[_lrr] + f[_rll] - f[_rlr] - f[_rrl] + f[_rrr];
    m[_212] = -f[_lll] - f[_llr] + f[_lrl] + f[_lrr] - f[_rll] - f[_rlr] + f[_rrl] + f[_rrr];
    m[_220] = f[_lll] + f[_llo] + f[_llr] + f[_lrl] + f[_lro] + f[_lrr] + f[_rll] + f[_rlo] + f[_rlr] + f[_rrl] +
              f[_rro] + f[_rrr];
    m[_221] = -f[_lll] + f[_llr] - f[_lrl] + f[_lrr] - f[_rll] + f[_rlr] - f[_rrl] + f[_rrr];
    m[_222] = f[_lll] + f[_llr] + f[_lrl] + f[_lrr] + f[_rll] + f[_rlr] + f[_rrl] + f[_rrr];
}

void raw_moment_to_pdf(double m[Q], double f[Q]) {
    f[_lll] = -0.125 * (m[_111] - m[_112] - m[_121] + m[_122] - m[_211] + m[_212] + m[_221] - m[_222]);
    f[_llo] = 0.25 * (m[_110] - m[_112] - m[_120] + m[_122] - m[_210] + m[_212] + m[_220] - m[_222]);
    f[_llr] = 0.125 * (m[_111] + m[_112] - m[_121] - m[_122] - m[_211] - m[_212] + m[_221] + m[_222]);
    f[_lol] = 0.25 * (m[_101] - m[_102] - m[_121] + m[_122] - m[_201] + m[_202] + m[_221] - m[_222]);
    f[_loo] = -0.5 * (m[_100] - m[_102] - m[_120] + m[_122] - m[_200] + m[_202] + m[_220] - m[_222]);
    f[_lor] = -0.25 * (m[_101] + m[_102] - m[_121] - m[_122] - m[_201] - m[_202] + m[_221] + m[_222]);
    f[_lrl] = 0.125 * (m[_111] - m[_112] + m[_121] - m[_122] - m[_211] + m[_212] - m[_221] + m[_222]);
    f[_lro] = -0.25 * (m[_110] - m[_112] + m[_120] - m[_122] - m[_210] + m[_212] - m[_220] + m[_222]);
    f[_lrr] = -0.125 * (m[_111] + m[_112] + m[_121] + m[_122] - m[_211] - m[_212] - m[_221] - m[_222]);
    f[_oll] = 0.25 * (m[_011] - m[_012] - m[_021] + m[_022] - m[_211] + m[_212] + m[_221] - m[_222]);
    f[_olo] = -0.5 * (m[_010] - m[_012] - m[_020] + m[_022] - m[_210] + m[_212] + m[_220] - m[_222]);
    f[_olr] = -0.25 * (m[_011] + m[_012] - m[_021] - m[_022] - m[_211] - m[_212] + m[_221] + m[_222]);
    f[_ool] = -0.5 * (m[_001] - m[_002] - m[_021] + m[_022] - m[_201] + m[_202] + m[_221] - m[_222]);
    f[_ooo] = 1.0 * (m[_000] - m[_002] - m[_020] + m[_022] - m[_200] + m[_202] + m[_220] - m[_222]);
    f[_oor] = 0.5 * (m[_001] + m[_002] - m[_021] - m[_022] - m[_201] - m[_202] + m[_221] + m[_222]);
    f[_orl] = -0.25 * (m[_011] - m[_012] + m[_021] - m[_022] - m[_211] + m[_212] - m[_221] + m[_222]);
    f[_oro] = 0.5 * (m[_010] - m[_012] + m[_020] - m[_022] - m[_210] + m[_212] - m[_220] + m[_222]);
    f[_orr] = 0.25 * (m[_011] + m[_012] + m[_021] + m[_022] - m[_211] - m[_212] - m[_221] - m[_222]);
    f[_rll] = 0.125 * (m[_111] - m[_112] - m[_121] + m[_122] + m[_211] - m[_212] - m[_221] + m[_222]);
    f[_rlo] = -0.25 * (m[_110] - m[_112] - m[_120] + m[_122] + m[_210] - m[_212] - m[_220] + m[_222]);
    f[_rlr] = -0.125 * (m[_111] + m[_112] - m[_121] - m[_122] + m[_211] + m[_212] - m[_221] - m[_222]);
    f[_rol] = -0.25 * (m[_101] - m[_102] - m[_121] + m[_122] + m[_201] - m[_202] - m[_221] + m[_222]);
    f[_roo] = 0.5 * (m[_100] - m[_102] - m[_120] + m[_122] + m[_200] - m[_202] - m[_220] + m[_222]);
    f[_ror] = 0.25 * (m[_101] + m[_102] - m[_121] - m[_122] + m[_201] + m[_202] - m[_221] - m[_222]);
    f[_rrl] = -0.125 * (m[_111] - m[_112] + m[_121] - m[_122] + m[_211] - m[_212] + m[_221] - m[_222]);
    f[_rro] = 0.25 * (m[_110] - m[_112] + m[_120] - m[_122] + m[_210] - m[_212] + m[_220] - m[_222]);
    f[_rrr] = 0.125 * (m[_111] + m[_112] + m[_121] + m[_122] + m[_211] + m[_212] + m[_221] + m[_222]);
}

void raw_moment_to_central_moment(double m[Q], double U[D], double k[Q]) {
    double u = U[0], v = U[1], w = U[2];
    double u2 = u * u, v2 = v * v, w2 = w * w;
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

void central_moment_to_raw_moment(double k[Q], double U[D], double m[Q]) {
    double u = U[0], v = U[1], w = U[2];
    double u2 = u * u, v2 = v * v, w2 = w * w;
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

void central_moment_to_cumulant(double k[Q], double U[D], double c[Q]) {
    double rho = k[_000];
    double u = U[0], v = U[1], w = U[2];
    c[_000] = rho;
    c[_100] = rho * u;
    c[_010] = rho * v;
    c[_001] = rho * w;
    c[_110] = k[_110];
    c[_101] = k[_101];
    c[_011] = k[_011];
    c[_200] = k[_200];
    c[_020] = k[_020];
    c[_002] = k[_002];
    c[_120] = k[_120];
    c[_102] = k[_102];
    c[_012] = k[_012];
    c[_210] = k[_210];
    c[_201] = k[_201];
    c[_021] = k[_021];
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
                         2 * (k[_110] * k[_102] + k[_011] * k[_201])) /
                            rho;
    c[_221] = k[_221] - (k[_200] * k[_021] + k[_020] * k[_201] + 4 * k[_110] * k[_111] +
                         2 * (k[_101] * k[_120] + k[_011] * k[_210])) /
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

void cumulant_to_central_moment(double c[Q], double shift[D], double k[Q]) {
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
    k[_120] = c[_120];
    k[_102] = c[_102];
    k[_012] = c[_012];
    k[_210] = c[_210];
    k[_201] = c[_201];
    k[_021] = c[_021];
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
                         2 * (k[_110] * k[_102] + k[_011] * k[_201])) /
                            rho;
    k[_221] = c[_221] + (k[_200] * k[_021] + k[_020] * k[_201] + 4 * k[_110] * k[_111] +
                         2 * (k[_101] * k[_120] + k[_011] * k[_210])) /
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

void pdf_to_cumulant(double f[Q], double c[Q]) {
    double m[Q], k[Q], U[D], rho;
    pdf_to_macroscopics(f, U, rho);
    pdf_to_raw_moment(f, m);
    raw_moment_to_central_moment(m, U, k);
    central_moment_to_cumulant(k, U, c);
}

void cumulant_to_pdf(double c[Q], double shift[D], double f[Q]) {
    double m[Q], k[Q], U[D];
    U[0] = c[_100] / c[_000];
    U[1] = c[_010] / c[_000];
    U[2] = c[_001] / c[_000];
    cumulant_to_central_moment(c, shift, k);
    central_moment_to_raw_moment(k, U, m);
    raw_moment_to_pdf(m, f);
}

void relax_cumulant(double c[Q], double omega, double c_post[Q]) {
    double rho = c[_000];
    c_post[_000] = c[_000];
    c_post[_100] = c[_100];
    c_post[_010] = c[_010];
    c_post[_001] = c[_001];
    c_post[_110] = (1 - omega) * c[_110];
    c_post[_101] = (1 - omega) * c[_101];
    c_post[_011] = (1 - omega) * c[_011];
    c_post[_200] = (rho + (1 - omega) * (2 * c[_200] - c[_020] - c[_002])) / 3;
    c_post[_020] = (rho + (1 - omega) * (-c[_200] + 2 * c[_020] - c[_002])) / 3;
    c_post[_002] = (rho + (1 - omega) * (-c[_200] - c[_020] + 2 * c[_002])) / 3;
    c_post[_120] = 0;
    c_post[_102] = 0;
    c_post[_012] = 0;
    c_post[_210] = 0;
    c_post[_201] = 0;
    c_post[_021] = 0;
    c_post[_111] = 0;
    c_post[_220] = 0;
    c_post[_202] = 0;
    c_post[_022] = 0;
    c_post[_211] = 0;
    c_post[_121] = 0;
    c_post[_112] = 0;
    c_post[_221] = 0;
    c_post[_212] = 0;
    c_post[_122] = 0;
    c_post[_222] = 0;
}

void compute_cumulant(double f[][Q], double c[][Q], int size[D]) {
    int len = size[0] * size[1] * size[2];
    for (int i = 0; i < len; i++) {
        pdf_to_cumulant(f[i], c[i]);
    }
}

void compute_collision(double c[][Q], double omega, double c_post[][Q], int size[D]) {
    int len = size[0] * size[1] * size[2];
    for (int i = 0; i < len; i++) {
        relax_cumulant(c[i], omega, c_post[i]);
    }
}

void compute_post_collision_pdf(double c[][Q], double shift[][D], double f[][Q], int size[D]) {
    int len = size[0] * size[1] * size[2];
    for (int i = 0; i < len; i++) {
        cumulant_to_pdf(c[i], shift[i], f[i]);
    }
}

void compute_streaming(double f_post[][Q], double f[][Q], int size[D]) {
    for (int i = 1; i < size[0] - 1; i++) {
        for (int j = 1; j < size[1] - 1; j++) {
            for (int k = 1; k < size[2] - 1; k++) {
                for (int q = 0; q < Q; q++) {
                    int src = index(i - Ve[q][0], j - Ve[q][1], k - Ve[q][2], size);
                    int dst = index(i, j, k, size);
                    f[dst][q] = f_post[src][q];
                }
            }
        }
    }
}

void apply_bc(double f[][Q], double f_post[][Q], double u_lid, int size[D]) {
    // left wall
    for (int j = 1; j < size[1] - 1; j++) {
        for (int k = 1; k < size[2] - 1; k++) {
            int id = index(1, j, k, size);
            int qlist[] = {_rll, _rlo, _rlr, _rol, _roo, _ror, _rrl, _rro, _rrr};
            for (int q : qlist) {
                f[id][q] = f_post[id][get_link_q(q)];
            }
        }
    }
    // right wall
    for (int j = 1; j < size[1] - 1; j++) {
        for (int k = 1; k < size[2] - 1; k++) {
            int id = index(size[0] - 2, j, k, size);
            int qlist[] = {_lll, _llo, _llr, _lol, _loo, _lor, _lrl, _lro, _lrr};
            for (int q : qlist) {
                f[id][q] = f_post[id][get_link_q(q)];
            }
        }
    }
    // bottom wall
    for (int i = 1; i < size[0] - 1; i++) {
        for (int j = 1; j < size[1] - 1; j++) {
            int id = index(i, j, 1, size);
            int qlist[] = {_llr, _lor, _lrr, _olr, _oor, _orr, _rlr, _ror, _rrr};
            for (int q : qlist) {
                f[id][q] = f_post[id][get_link_q(q)];
            }
        }
    }
    // top moving wall
    for (int i = 1; i < size[0] - 1; i++) {
        for (int j = 1; j < size[1] - 1; j++) {
            int id = index(i, j, size[2] - 2, size);
            int qlist[] = {_lll, _lol, _lrl, _oll, _ool, _orl, _rll, _rol, _rrl};
            for (int q : qlist) {
                int link = get_link_q(q);
                f[id][q] = f_post[id][link] - Ve[link][0] * 2 * u_lid * Wght[link] * CSQI;
            }
        }
    }
    // back wall
    for (int i = 1; i < size[0] - 1; i++) {
        for (int k = 1; k < size[2] - 1; k++) {
            int id = index(i, 1, k, size);
            int qlist[] = {_lrl, _lro, _lrr, _orl, _oro, _orr, _rrl, _rro, _rrr};
            for (int q : qlist) {
                f[id][q] = f_post[id][get_link_q(q)];
            }
        }
    }
    // front wall
    for (int i = 1; i < size[0] - 1; i++) {
        for (int k = 1; k < size[2] - 1; k++) {
            int id = index(i, size[1] - 2, k, size);
            int qlist[] = {_lll, _llo, _llr, _oll, _olo, _olr, _rll, _rlo, _rlr};
            for (int q : qlist) {
                f[id][q] = f_post[id][get_link_q(q)];
            }
        }
    }
}

template <typename T>
void copy_array(T *src, T *dst, int len) {
    for (int i = 0; i < len; i++) {
        dst[i] = src[i];
    }
}

template <typename T, int N>
void copy_array(T src[][N], T dst[][N], int len) {
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < N; j++) {
            dst[i][j] = src[i][j];
        }
    }
}

class Cumu {
public:
    double (*f)[Q], (*f_post)[Q], (*c)[Q], (*c_post)[Q], (*shift)[D];
    int size[D];
    double tau;

    Cumu(int size[D], double tau) {
        copy_array(size, this->size, 3);
        this->tau = tau;

        int len = size[0] * size[1] * size[2];
        f = new double[len][Q]();
        f_post = new double[len][Q]();
        c = new double[len][Q]();
        c_post = new double[len][Q]();
    }

    ~Cumu() {
        delete[] f;
        delete[] f_post;
        delete[] c;
        delete[] c_post;
    }
};

Cumu *init(int size[D], double tau) {
    Cumu *cumu = new Cumu(size, tau);
    int len = size[0]*size[1]*size[2];
    for (int i = 0; i < len; i ++) {
        double shift[] = {0, 0, 0};
        double rho = 1;
        double U[] = {0, 0, 0};
        double ceq[Q];
        get_eq_cumulant(rho, U, ceq);
        copy_array(shift, cumu->shift[i], Q);
        cumulant_to_pdf(ceq, shift, cumu->f[i]);
    }

    return cumu;
}

void finalize(Cumu *cumu) {
    delete cumu;
}

