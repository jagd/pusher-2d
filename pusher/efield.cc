#include "efield.h"
#include <limits>
#include <cmath>

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
