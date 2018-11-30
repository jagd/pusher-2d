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


double LeapFrogPusher2D::pTheta(double z, double r, double uTheta) const
{
    return r*(uTheta*M0 + magfield_->aTheta(z, r)*Q0);
}


void LeapFrogPusher2D::step(double dt)
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
        Q0 / M0 * (ez - vTheta * magfield_->br(z, r))
        ,
        Q0 / M0 * er + vTheta*(uTheta / r + Q0 / M0 * magfield_->bz(z, r))
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


void LeapFrogPusher2D::setElectronInfo(
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

double LeapFrogPusher2D::gammaCurrent() const
{
    const double Ekin = totalEnergy_-Q0*efield_->pot(pos_.z, pos_.r);
//    assert(Ekin > = 0);
    const double gamma = 1.0 + Ekin/(M0*C0*C0);
    return gamma;
}

PV2D LeapFrogPusher2D::vLastHalf() const
{
    return uLastHalf_ / gammaCurrent();
}

PV2D LeapFrogPusher2D::pos() const
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
    const auto u1 = u_;
    const auto v1 = u1 / usqr2gamma(norm2(u1));
    const auto hkr1 = dt * v1;
	const auto hku1 = dt * a3D(pos_, v1);

    const auto u2 = u1 + hku1 / 2;
    const auto v2 = u2 / usqr2gamma(norm2(u2));
    const auto hkr2 = dt * v2;
	const auto hku2 = dt * a3D(pos_ + hkr1/2, v2);

    const auto u3 = u1 + hku2 / 2;
    const auto v3 = u3 / usqr2gamma(norm2(u3));
	const auto hkr3 = dt * v3;
	const auto hku3 = dt * a3D(pos_ + hkr2/2, v3);

    const auto u4 = u1 + hku3;
    const auto v4 = u4 / usqr2gamma(norm2(u4));
    const auto hkr4 = dt * v4;
	const auto hku4 = dt * a3D(pos_ + hkr3, v4);

    u_ = u1 + 1.0 / 6 * (hku1 + 2 * hku2 + 2 * hku3 + hku4);
	pos_ += 1.0 / 6 * (hkr1 + 2 * hkr2 + 2 * hkr3 + hkr4);
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

void LeapFrogPusher2DSync::step(double dt)
{
    const double dtHalf = dt / 2;
    // the first half step: use a correct v
    const auto v = u_ / gammaCurrent();
    const auto posInter = pos_ + dtHalf * v;
    const double p2mrInter = pTheta_/ (M0*posInter.r);
    const double uThetaInter =  p2mrInter - Q0 / M0 * magfield_->aTheta(posInter.z, posInter.r);
    const double gInter = gammaAt(posInter);
    const double vThetaInter = uThetaInter / gInter;
    const double ezInter = efield_->ez(posInter.z, posInter.r);
    const double erInter = efield_->er(posInter.z, posInter.r);
    const PV2D dudtInter(
        Q0 / M0 * (ezInter - vThetaInter * magfield_->br(posInter.z, posInter.r))
        ,
        Q0 / M0 * erInter +
        vThetaInter*(uThetaInter / posInter.r + Q0 / M0 * magfield_->bz(posInter.z, posInter.r))
    );
    const auto uNext = u_ + dudtInter * dt;
    // discriminant
    const double disc = gInter * gInter + (2 * Q0 / M0 / C0 / C0) * dt * (
        ezInter*uNext.z + erInter*uNext.r
    );
    const double gNextAppr = (gInter + std::sqrt(disc)) / 2;
    const auto vNextAppr = uNext / gNextAppr;
    // the second half step: use a approximated (extrapolated) v
    pos_ = posInter + vNextAppr * dtHalf;
}

void LeapFrogPusher2DSync::setElectronInfo(double z, double r, double uz, double ur, double uTheta)
{
    pos_ = PV2D(z, r);
    u_ = PV2D(uz, ur);
    pTheta_ = pTheta(z, r, uTheta);
    const double gamma = std::sqrt(1 + (uz*uz + ur * ur + uTheta * uTheta) / (C0*C0));
    totalEnergy_ = (gamma - 1) * (M0*C0*C0) + Q0*efield_->pot(z, r);
}

double LeapFrogPusher2DSync::pTheta(double z, double r, double uTheta) const
{
    return r*(uTheta*M0 + magfield_->aTheta(z, r)*Q0);
}

double LeapFrogPusher2DSync::gammaAt(const PV2D & pos) const
{
    const double Ekin = totalEnergy_-Q0*efield_->pot(pos.z, pos.r);
//    assert(Ekin > = 0);
    const double gamma = 1.0 + Ekin/(M0*C0*C0);
    return gamma;
}
