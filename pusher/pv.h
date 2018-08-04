#ifndef WU_PV_H__
#define WU_PV_H__

#include <iosfwd>
#include <cmath>

/**@breif Point or Vector in 2D cylindrical coordinate */
struct PV2D {
    PV2D() = default;
    PV2D(double z, double r);
	double z;
	double r;
};

/**@breif Point or Vector in 3D Cartesian coordinate */
struct PV3D {
    PV3D() = default;
    PV3D(double x, double y, double z) : x(x), y(y), z(z) {}
	double x;
	double y;
	double z;
};

// inline implementations
inline PV2D fromPV3D(const PV3D &p)
{
    return PV2D(p.z, std::sqrt(p.x*p.x+p.y*p.y));
}

inline double dot(const PV2D& a, const PV2D& b)
{
    return a.r*b.r + a.z*b.z;
}

inline PV2D &operator+=(PV2D &a, const PV2D &b)
{
    a.z += b.z;
    a.r += b.r;
    return a;
}

inline PV2D &operator*=(PV2D &v, double a)
{
    v.z *= a;
    v.r *= a;
    return v;
}

inline PV2D &operator/=(PV2D &v, double a)
{
    v.z /= a;
    v.r /= a;
    return v;
}

inline PV2D operator*(double a, const PV2D &v)
{
    auto x = v;
    x *= a;
    return x;
}

inline PV2D operator*(const PV2D &v, double a)
{
    return a*v;
}

inline PV2D operator/(const PV2D &v, double a)
{
    auto x = v;
    x /= a;
    return x;
}


inline PV3D cross(const PV3D& a, const PV3D& b)
{
    return PV3D(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}

inline double dot(const PV3D& a, const PV3D& b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline PV3D &operator/=(PV3D &v, double a)
{
    v.x /= a;
    v.y /= a;
    v.z /= a;
    return v;
}

inline PV3D operator/(const PV3D &v, double a)
{
    auto x = v;
    x /= a;
    return x;
}

inline PV3D &operator*=(PV3D &v, double a)
{
    v.x *= a;
    v.y *= a;
    v.z *= a;
    return v;
}

inline PV3D &operator+=(PV3D &a, const PV3D &b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

inline PV3D &operator-=(PV3D &a, const PV3D &b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

inline PV3D operator*(double a, const PV3D &v)
{
    auto x = v;
    x *= a;
    return x;
}

inline PV3D operator*(const PV3D &v, double a)
{
    return a*v;
}

inline PV3D operator+(const PV3D &a, const PV3D &b)
{
    auto x = a;
    x += b;
    return x;
}

inline PV3D operator-(const PV3D &a, const PV3D &b)
{
    auto x = a;
    x -= b;
    return x;
}

std::ostream &operator<<(std::ostream &, const PV2D &);
std::ostream &operator<<(std::ostream &, const PV3D &);

#endif
