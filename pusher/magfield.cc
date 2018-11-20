#include "magfield.h"
#include <cmath>
#include <cassert>

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

BiUniformMagField::BiUniformMagField(double z0, double b1, double b2)
    : z0_(z0), z02_(z0_*z0), z03_(z02_*z0), z04_(z03_*z0), z05_(z04_*z0),
      b1_(b1), b2_(b2), coef_(15.0/16.0*(b2-b1)/z0)
{
    assert(z0 >= 0);
}

double BiUniformMagField::aTheta(double z, double r) const
{
    return bz(z, r)*r/2;
}

double BiUniformMagField::aTheta2z(double z, double r) const
{
    if (std::abs(z) > z0_) {
        return 0;
    }

    return coef_ * (std::pow(z, 4) / z04_ - 2 * std::pow(z, 2) / z02_ + 1)*r/2;
}

double BiUniformMagField::aTheta2r(double z, double r) const
{
    return bz(z, r)/2;
}

double BiUniformMagField::br(double z, double r) const
{
    if (std::abs(z) > z0_)
        return 0.0;

    const double nz = z / z0_;
    const double a = 1 - nz * nz;
    return -0.5*coef_*r*(a*a);
}

double BiUniformMagField::bz(double z, double /*r*/) const
{
    if (z <= -z0_) {
        return b1_;
    } else if (z >= z0_) {
        return b2_;
    }

    return coef_ * (1.0 / 5 * (std::pow(z, 5) + z05_) / z04_ - 2.0 / 3 * (std::pow(z, 3) + z03_) / z02_ + z + z0_) + b1_;
}
