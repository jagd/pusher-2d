#include "efield.h"
#include <limits>
#include <cmath>
#include <cassert>

double IStaticEField::pot(double, double) const
{
    return std::numeric_limits<double>::quiet_NaN();
}


IStaticEField::~IStaticEField()
{
}


MirroredEzField::MirroredEzField(double ez): ez_(ez)
{
}


double MirroredEzField::pot(double z, double) const
{
    return -ez_ * std::abs(z);
}


double MirroredEzField::ez(double z, double) const
{
    return (z >= 0) ? ez_ : -ez_;
}


double MirroredEzField::er(double, double) const
{
    return 0;
}

ConstEzField::ConstEzField(double ez): ez_(ez)
{
}

double ConstEzField::pot(double z, double) const
{
    return -z*ez_;
}

double ConstEzField::ez(double, double) const
{
    return ez_;
}

double ConstEzField::er(double, double) const
{
    return 0;
}

ConstErField::ConstErField(double er) : er_(er)
{
}

double ConstErField::pot(double, double r) const
{
    assert(r >= 0);
    return -r*er_;
}

double ConstErField::ez(double, double) const
{
    return 0;
}

double ConstErField::er(double, double) const
{
    return er_;
}

LinearEzField::LinearEzField(double ez0, double k): ez0_(ez0), k_(k)
{
}

double LinearEzField::pot(double z, double) const
{
    return -(ez0_ * z + k_ * z*z / 2);
}

double LinearEzField::ez(double z, double) const
{
    return ez0_ + k_ * z;
}

double LinearEzField::er(double, double) const
{
    return 0.0;
}
