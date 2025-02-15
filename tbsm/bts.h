#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <cstdio>

class BtsReader {
public:
    const int nv = 3;
    int16_t magic;
    int32_t ny, nz, ntwr, ntime;
    float dz, dy, dt, uhub, zhub, zlow;
    float scl[3], off[3];
    int32_t nchar;
    char *info;
    int16_t *buffer;

    std::ifstream ifs;
    size_t header;

    ~BtsReader() {
        delete[] info;
        delete[] buffer;
        ifs.close();
    }

    int btsid(int j, int k, int l) {
        return k*ny*nv + j*nv + l;
    }

    BtsReader(std::string path) {
        ifs.open(path, std::ios::binary);
        ifs.read((char*)&magic, sizeof(int16_t));
        ifs.read((char*)&nz, sizeof(int32_t));
        ifs.read((char*)&ny, sizeof(int32_t));
        ifs.read((char*)&ntwr, sizeof(int32_t));
        ifs.read((char*)&ntime, sizeof(int32_t));
        ifs.read((char*)&dz, sizeof(float));
        ifs.read((char*)&dy, sizeof(float));
        ifs.read((char*)&dt, sizeof(float));
        ifs.read((char*)&uhub, sizeof(float));
        ifs.read((char*)&zhub, sizeof(float));
        ifs.read((char*)&zlow, sizeof(float));
        ifs.read((char*)&scl[0], sizeof(float));
        ifs.read((char*)&off[0], sizeof(float));
        ifs.read((char*)&scl[1], sizeof(float));
        ifs.read((char*)&off[1], sizeof(float));
        ifs.read((char*)&scl[2], sizeof(float));
        ifs.read((char*)&off[2], sizeof(float));
        ifs.read((char*)&nchar, sizeof(int32_t));
        info = new char[nchar];
        ifs.read(info, sizeof(char)*nchar);
        buffer = new int16_t[ny*nz*nv];

        header = sizeof(int16_t) + 5*sizeof(int32_t) + 12*sizeof(float) + nchar*sizeof(char);
    }

    void print_info() {
        printf("BTS FILE INFO\n");
        printf("\tid = %d\n", magic);
        printf("\tsize = (%d %d)\n", ny, nz);
        printf("\tgrid size = (%f %f)\n", dy, dz);
        printf("\tdt = %f\n", dt);
        printf("\ttower = %d\n", ntwr);
        printf("\tstep = %d\n", ntime);
        printf("\tu mean at hub height = %f\n", uhub);
        printf("\thub height = %f\n", zhub);
        printf("\tbottom height = %f\n", zlow);
        printf("\tscale = (%f %f %f)\n", scl[0], scl[1], scl[2]);
        printf("\toffset = (%f %f %f)\n", off[0], off[1], off[2]);
        printf("\t%s\n", info);
    }

    template<typename T>
    T peek(int t, int j, int k, int l) {
        size_t skip = header + t*nv*(ny*nz + ntwr)*sizeof(int16_t) + btsid(j, k, l)*sizeof(int16_t);
        int16_t value;
        ifs.seekg(skip);
        ifs.read((char*)&value, sizeof(int16_t));
        return ((T)value - off[l])/scl[l];
    }

    template<typename T>
    void load_to_array(T *ptr, int t) {
        if (t >= ntime) {
            printf("BTS READER ERROR: TIME OUT OF BOUND! %d > %d\n", t, ntime - 1);
        }
        size_t skip = header + t*nv*(ny*nz + ntwr)*sizeof(int16_t);
        size_t len  = nv*ny*nz*sizeof(int16_t);
        ifs.seekg(skip, std::ios::beg);
        ifs.read((char*)buffer, len);
        for (int k = 0; k < nz; k ++) {
        for (int j = 0; j < ny; j ++) {
        for (int l = 0; l < nv; l ++) {
            ptr[btsid(j, k, l)] = (buffer[btsid(j, k, l)] - off[l])/scl[l];
        }}}
    }
};