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

////////////////////////////////////////////////////////////////////////////////

double LinearBzField::aTheta(double z, double r) const
{
	return bz(z, r) * r / 2;
}

double LinearBzField::aTheta2z(double, double r) const
{
	return k_ * r / 2;
}

double LinearBzField::aTheta2r(double z, double r) const
{
	return bz(z, r) / 2;
}

double LinearBzField::br(double, double r) const
{
	return -r * (k_ / 2);
}

double LinearBzField::bz(double z, double) const
{
	return b0_ + k_ * z;
}
