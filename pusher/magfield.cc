#include "magfield.h"

IStaticMagField::IStaticMagField()
{

}

IStaticMagField::~IStaticMagField()
{}

double ConstBzField::br(double, double) const
{
    return 0;
}

double ConstBzField::bz(double, double) const
{
    return bz_;
}

double ConstBzField::aTheta(double, double r) const
{
    return bz_*r/2;
}

double ConstBzField::aTheta2z(double, double) const
{
    return 0;
}

double ConstBzField::aTheta2r(double, double) const
{
    return bz_/2;
}
