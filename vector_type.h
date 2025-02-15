#pragma once

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

    vector_t<T, N> operator-() {
        vector_t<T, N> tmp;
        for (int i = 0; i < N; i ++) {
            tmp[i] = - m[i];
        } 
        return tmp;
    }
};

template<typename T, int N>
vector_t<T, N> operator*(const T &s, const vector_t<T, N> &v) {
    vector_t<T, N> tmp;
    for (int i = 0; i < N; i ++) {
        tmp[i] = s*v[i];
    }
    return tmp;
}