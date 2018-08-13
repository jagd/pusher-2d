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


MirOscZEField::MirOscZEField(double ez): ez_(ez)
{
}


double MirOscZEField::pot(double z, double) const
{
    return -ez_ * std::abs(z);
}


double MirOscZEField::ez(double z, double) const
{
    return (z >= 0) ? ez_ : -ez_;
}


double MirOscZEField::er(double, double) const
{
    return 0;
}
