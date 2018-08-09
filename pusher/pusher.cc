#include "pusher.h"
#include "pv.h"
#include <cmath>

IPusher2D::IPusher2D(
    std::shared_ptr<IStaticEField> e,
    std::shared_ptr<IStaticMagField> m
): efield_(e), magfield_(m)
{
}


IPusher2D::~IPusher2D()
{
}


void BorisPusher::setElectronInfo(
    double x, double y, double z,
    double ux, double uy, double uz
)
{
    pos_ = PV3D(x, y, z);
    uLastHalf_ = PV3D(ux, uy, uz);
}


void BorisPusher::step(double dt)
{
    const PV2D zr = fromPV3D(pos_);
    const PV3D posNorm = (zr.r == 0) ? PV3D(0,0,1) : (pos_/std::sqrt(dot(pos_, pos_)));
    const auto er = efield_->er(zr.z, zr.r);
    const PV3D e3d(er*posNorm.x, er*posNorm.y, efield_->ez(zr.z, zr.r));
    const PV3D uMinus = uLastHalf_ + (Q0/2*M0)*dt*e3d;
    const double uMinusNorm2 = dot(uMinus, uMinus);
    const double gammaN = std::sqrt(1.0 + uMinusNorm2/(C0*C0));

    const auto br = magfield_->br(zr.z, zr.r);
    const PV3D b3d(br*posNorm.x, br*posNorm.y, magfield_->bz(zr.z, zr.r));
    const PV3D t = (Q0/M0/2)*dt/gammaN*b3d;
    const PV3D s = 2.0/(1+dot(t,t)) * t;
    const PV3D uPrime = uMinus + cross(uMinus, t);
    const PV3D uPlus = uMinus + cross(uPrime, s);

    const PV3D uNextHalf = uPlus + (Q0/2*M0)*dt*e3d;
    pos_ = pos_ + dt/(std::sqrt(1+dot(uNextHalf, uNextHalf)/(C0*C0))) * uNextHalf;
    uLastHalf_ = uNextHalf;
}


PV3D BorisPusher::pos() const
{
    return pos_;
}


PV3D BorisPusher::u() const
{
    return uLastHalf_;
}


PV3D BorisPusher::gamma() const
{
    return gammaAtPos_;
}

double v2gamma(double v)
{
    return 1.0/(std::sqrt(1.0-v*v/(C0*C0)));
}

double v2u(double v)
{
    return v2gamma(v)*v;
}

double u2gamma(double u)
{
    return std::sqrt(1.0+u*u/(C0*C0));
}
