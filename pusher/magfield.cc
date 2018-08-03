#include "magfield.h"

IStaticMagField::IStaticMagField()
{

}

IStaticMagField::~IStaticMagField()
{}

double HomogeneousMagField::br(double, double) const
{
    return 0;
}

double HomogeneousMagField::bz(double, double) const
{
    return bz_;
}

double HomogeneousMagField::aphi(double, double r) const
{
    return bz_*r/2;
}

double HomogeneousMagField::aphi2z(double, double) const
{
    return 0;
}

double HomogeneousMagField::aphi2r(double, double) const
{
    return bz_/2;
}
