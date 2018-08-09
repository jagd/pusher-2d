#include "efield.h"
#include <limits>

IStaticEField::IStaticEField()
{

}

double IStaticEField::pot(double, double) const
{
    return std::numeric_limits<double>::quiet_NaN();
}

IStaticEField::~IStaticEField()
{
}
