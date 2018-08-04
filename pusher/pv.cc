#include "pv.h"
#include <ostream>

PV2D::PV2D(double z, double r): z(z), r(r)
{
	//  assert(r >= 0); is ignored since vectors may have negative r
}

std::ostream &operator<<(std::ostream &os, const PV2D &c)
{
    os << c.z << ' ' << c.r;
    return os;
}

std::ostream &operator<<(std::ostream &os, const PV3D &c)
{
    os << c.x << ' ' << c.y << ' ' << c.z;
    return os;
}

