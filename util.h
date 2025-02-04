template<typename T>
T bilinear_interpolate(
    T v0, T v1, T v2, T v3,
    T x0, T x1, T y0, T y1,
    T xx, T yy
) {
    T ax = (xx - x0)/(x1 - x0);
    T ay = (yy - y0)/(y1 - y0);
    T v4 = (1 - ax)*v0 + ax*v1;
    T v5 = (1 - ax)*v2 + ax*v3;
    T vc = (1 - ay)*v4 + ay*v5;
    return vc;
}
