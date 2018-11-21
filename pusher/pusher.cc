#include "pusher.h"
#include "pv.h"
#include <limits>
#include <cmath>
#include <cassert>

double gamma2v(double g)
{
    return C0 * std::sqrt(1-1/(g*g));
}

double v2gamma(double v)
{
    return 1.0/(std::sqrt(1.0-v*v/(C0*C0)));
}

double vsqr2gamma(double vsqr)
{
    return 1.0/(std::sqrt(1.0-vsqr/(C0*C0)));
}

double v2u(double v)
{
    return v2gamma(v)*v;
}


double u2gamma(double u)
{
    return std::sqrt(1.0+u*u/(C0*C0));
}

double usqr2gamma(double usqr)
{
    return std::sqrt(1.0+usqr/(C0*C0));
}


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
	gammaAtPos_ = std::numeric_limits<double>::quiet_NaN();
}


void BorisPusher::step(double dt)
{
    const PV2D zr = fromPV3D(pos_);
    const PV3D planarNorm = (zr.r == 0) ? PV3D(0,0,1) : (pos_/zr.r);
    const auto er = efield_->er(zr.z, zr.r);
    const PV3D e3d(er*planarNorm.x, er*planarNorm.y, efield_->ez(zr.z, zr.r));
    const PV3D uMinus = uLastHalf_ + (Q0/2/M0)*dt*e3d;
    const double uMinusNorm2 = norm2(uMinus);
    const double gammaN = std::sqrt(1.0 + uMinusNorm2/(C0*C0));

    const auto br = magfield_->br(zr.z, zr.r);
    const PV3D b3d(br*planarNorm.x, br*planarNorm.y, magfield_->bz(zr.z, zr.r));
    const PV3D t = (Q0/M0/2)*dt/gammaN*b3d;
    const PV3D s = 2.0/(1+norm2(t)) * t;
    const PV3D uPrime = uMinus + cross(uMinus, t);
    const PV3D uPlus = uMinus + cross(uPrime, s);

    const PV3D uNextHalf = uPlus + (Q0/2/M0)*dt*e3d;
    pos_ = pos_ + dt/(std::sqrt(1+norm2(uNextHalf)/(C0*C0))) * uNextHalf;
    uLastHalf_ = uNextHalf;
    gammaAtPos_ = gammaN;
}


PV3D BorisPusher::pos() const
{
    return pos_;
}


PV3D BorisPusher::uLastHalf() const
{
    return uLastHalf_;
}


double BorisPusher::gammaCurrent() const
{
    return gammaAtPos_;
}


double APhiPusher::pTheta(double z, double r, double uTheta) const
{
    return r*(uTheta*M0 + magfield_->aTheta(z, r)*Q0);
}


void APhiPusher::step(double dt)
{
    const double g = gammaCurrent();
    const double z = pos_.z;
    const double r = pos_.r;
    const double ez = efield_->ez(z, r);
    const double er = efield_->er(z, r);
    const double p2mr = pTheta_/ (M0*r);
    const double uTheta = p2mr - Q0 / M0 * magfield_->aTheta(z, r);
    const double vTheta = uTheta / g;
    const PV2D dudt(
        Q0 / M0 * (vTheta* magfield_->aTheta2z(z, r) + ez)
        ,
        vTheta*(p2mr / r + Q0 / M0 * magfield_->aTheta2r(z, r)) + Q0 / M0 * er
    );
    const PV2D uNextHalf = uLastHalf_ + dt * dudt;
#ifdef GAMMA_CORRECTION_APHI
    // discriminant
    const double disc = g * g + (2 * Q0 / M0 / C0 / C0) * dt * (
        ez*uNextHalf.z + er*uNextHalf.r
    );
    const double gammaNextHalf = (g + std::sqrt(disc)) / 2;
    pos_ += uNextHalf / gammaNextHalf * dt;
#else // no GAMMA_CORRECTION_APHI
    pos_ += uNextHalf / g * dt;
#endif
    uLastHalf_ = uNextHalf;
}


void APhiPusher::setElectronInfo(
    double z,
    double r,
    double uzLastHalf,
    double urLastHalf,
    double pTheta,
    double gamma
)
{
    assert(gamma >= 1);
    pos_ = PV2D(z, r);
    uLastHalf_ = PV2D(uzLastHalf, urLastHalf);
    pTheta_ = pTheta;
    totalEnergy_ = (gamma - 1) * (M0*C0*C0) + Q0*efield_->pot(z, r);
}

double APhiPusher::gammaCurrent() const
{
    const double Ekin = totalEnergy_-Q0*efield_->pot(pos_.z, pos_.r);
//    assert(Ekin > = 0);
    const double gamma = 1.0 + Ekin/(M0*C0*C0);
    return gamma;
}

PV2D APhiPusher::vLastHalf() const
{
    return uLastHalf_ / gammaCurrent();
}

PV2D APhiPusher::pos() const
{
    return pos_;
}

void LeapFrogPusher::setElectronInfo(double x, double y, double z, double ux, double uy, double uz)
{
	pos_ = PV3D(x, y, z);
	uLastHalf_ = PV3D(ux, uy, uz);
}

void LeapFrogPusher::step(double dt)
{
    const PV2D zr = fromPV3D(pos_);
    const PV3D planarNorm = (zr.r == 0) ? PV3D(0,0,1) : (pos_/zr.r);
    const auto br = magfield_->br(zr.z, zr.r);
    const PV3D b3d(br*planarNorm.x, br*planarNorm.y, magfield_->bz(zr.z, zr.r));
    const auto er = efield_->er(zr.z, zr.r);
    const PV3D e3d(er*planarNorm.x, er*planarNorm.y, efield_->ez(zr.z, zr.r));
    const double gammaLastHalf = usqr2gamma(norm2(uLastHalf_));
    const PV3D vLastHalf = uLastHalf_ / gammaLastHalf;
    const PV3D a = Q0/M0 * (cross(vLastHalf, b3d) + e3d);
    const PV3D uNextHalf = uLastHalf_ + a * dt; // becomes uNextHalf
    const PV3D vNextHalf = uNextHalf / usqr2gamma(norm2(uNextHalf));
	pos_ += vNextHalf * dt;
	uLastHalf_ = uNextHalf;
}


PV3D LeapFrogPusher::pos() const
{
	return pos_;
}

PV3D LeapFrogPusher::uLastHalf() const
{
	return uLastHalf_;
}

double LeapFrogPusher::gammaLastHalf() const
{
    return usqr2gamma(norm2(uLastHalf_));
}

void RK4Pusher::setElectronInfo(double x, double y, double z, double ux, double uy, double uz)
{
	pos_ = PV3D(x, y, z);
	u_ = PV3D(ux, uy, uz);
}

void RK4Pusher::step(double dt)
{
    const auto u = u_;

    const auto vInit = u / usqr2gamma(norm2(u));
    const auto r0 = dt * vInit;
	const auto u0 = dt * a3D(pos_, vInit);

    const auto u0HalfMore = u + u0 / 2;
    const auto v0HalfMore = u0HalfMore / usqr2gamma(norm2(u0HalfMore));
    const auto r1 = dt * v0HalfMore;
	const auto u1 = dt * a3D(pos_ + r0/2, v0HalfMore);

    const auto u1HalfMore = u + u1 / 2;
    const auto v1HalfMore = u1HalfMore / usqr2gamma(norm2(u1HalfMore));
	const auto r2 = dt * v1HalfMore;
	const auto u2 = dt * a3D(pos_ + r1/2, v1HalfMore);

    const auto u2More = u + u2;
    const auto v2More = u2More / usqr2gamma(norm2(u2More));
    const auto r3 = dt * v2More;
	const auto u3 = dt * a3D(pos_ + r2, v2More);

    u_ = u + 1.0 / 6 * (u0 + 2 * u1 + 2 * u2 + u3);
	pos_ += 1.0 / 6 * (r0 + 2 * r1 + 2 * r2 + r3);
}

PV3D RK4Pusher::pos() const
{
	return pos_;
}

PV3D RK4Pusher::uCurrent() const
{
	return u_;
}

PV3D RK4Pusher::b3D(const PV3D & pos) const
{
    const PV2D zr = fromPV3D(pos);
    const PV3D planarNorm = (zr.r == 0) ? PV3D(0,0,1) : (pos/zr.r);
    const auto br = magfield_->br(zr.z, zr.r);
    return PV3D(br*planarNorm.x, br*planarNorm.y, magfield_->bz(zr.z, zr.r));
}

PV3D RK4Pusher::e3D(const PV3D & pos) const
{
    const PV2D zr = fromPV3D(pos);
    const PV3D planarNorm = (zr.r == 0) ? PV3D(0,0,1) : (pos/zr.r);
	const auto er = efield_->er(zr.z, zr.r);
	return PV3D(er*planarNorm.x, er*planarNorm.y, efield_->ez(zr.z, zr.r));
}

PV3D RK4Pusher::a3D(const PV3D & pos, const PV3D & v) const
{
    return Q0/M0 * (cross(v, b3D(pos)) + e3D(pos));
}

double RK4Pusher::gammaCurrent() const
{
    return usqr2gamma(norm2(u_));
}
