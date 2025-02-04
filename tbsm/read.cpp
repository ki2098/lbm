#include <cmath>
#include "bts.h"
#include "../vector_type.h"
#include "../util.h"

using namespace std;

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
    }

    ~InflowHandler() {
        delete[] inflow;
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
    }
};

int main() {
    int ny = 100, nz = 1;
    InflowHandler<2> inflow(ny, nz, "TurbSim.bts");
    inflow.print_info();
    float oy = inflow.oy, oz = inflow.oz;
    inflow.u_scale = inflow.uhub;
    
    inflow.apply_inflow(0.);
    float dy = inflow.dy, dz = inflow.dz;

    FILE *file = fopen("data/inflow.csv", "w");
    fprintf(file, "x,y,z,u,v,w\n");
    for (int k = 0; k < nz; k ++) {
    for (int j = 0; j < ny; j ++) {
        fprintf(
            file,
            "%f,%f,%f,%f,%f,%f\n",
            0., dy*j + oy, dz*k + oz,
            inflow.inflow[j*nz + k][0], inflow.inflow[j*nz + k][1], 0.
        );
    }}
    fclose(file);
}