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
#if 0
    const PV2D dudt(
        Q0 / M0 * (vTheta* magfield_->aTheta2z(z, r) + ez)
        ,
        vTheta*(p2mr / r + Q0 / M0 * magfield_->aTheta2r(z, r)) + Q0 / M0 * er
    );
#else
    const PV2D dudt(
        -vTheta * magfield_->br(z, r) + Q0 / M0 * ez
        ,
        vTheta*(uTheta / r + Q0 / M0 * magfield_->bz(z, r)) + Q0 / M0 * er
    );
#endif
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

/// Based on the current position, advance dt/2
double RK4Pusher2D::extGammaHalf(
    double baseGamma, //! avoid duplicated gamma calculation
    double dt,
    const PV2D &uTarget
)
{
    const auto z = pos_.z;
    const auto r = pos_.r;
    const double disc = baseGamma * baseGamma + (2 * Q0 / M0 / C0 / C0) * dt * (
        efield_->ez(z, r)*uTarget.r + efield_->er(z, r)*uTarget.r
    );
    return (baseGamma + std::sqrt(disc)) / 2;
}

void RK4Pusher2D::step(double dt)
{

    const double g1 = gammaAt(pos_);
    const auto u1 = u_;
    const auto v1 = u1 / g1;
    const auto hkr1 = dt * v1;
	const auto hku1 = dt * dudtAt(pos_);

    const auto u2 = u1 + hku1 / 2;
    const auto v2 = u2 / extGammaHalf(g1, dt, u2);
    const auto hkr2 = dt * v2;
	const auto hku2 = dt * dudtAt(pos_ + hkr1/2);

    const auto u3 = u1 + hku2 / 2;
    const auto v3 = u3 / extGammaHalf(g1, dt, u3);
	const auto hkr3 = dt * v3;
	const auto hku3 = dt * dudtAt(pos_ + hkr2/2);

    const auto u4 = u1 + hku3;
    const auto v4 = u4 / extGammaHalf(g1, 2*dt, u4);
    const auto hkr4 = dt * v4;
	const auto hku4 = dt * dudtAt(pos_ + hkr3);

    u_ = u1 + 1.0 / 6 * (hku1 + 2 * hku2 + 2 * hku3 + hku4);
	pos_ += 1.0 / 6 * (hkr1 + 2 * hkr2 + 2 * hkr3 + hkr4);
}

void RK4Pusher2D::setElectronInfo(double z, double r, double uz, double ur, double uTheta)
{
    pos_ = PV2D(z, r);
    u_ = PV2D(uz, ur);
    pTheta_ = r*(uTheta*M0 + magfield_->aTheta(z, r)*Q0);
    const double gamma = std::sqrt(
        1 + (uz*uz + ur * ur + uTheta * uTheta) / (C0*C0)
    );
    totalEnergy_ = (gamma - 1) * (M0*C0*C0) + Q0*efield_->pot(z, r);
}

PV2D RK4Pusher2D::pos() const
{
    return pos_;
}

double RK4Pusher2D::gammaCurrent() const
{
    return gammaAt(pos_);
}

double RK4Pusher2D::gammaAt(const PV2D &pos) const
{
    const double Ekin = totalEnergy_-Q0*efield_->pot(pos_.z, pos_.r);
//    assert(Ekin > = 0);
    const double gamma = 1.0 + Ekin/(M0*C0*C0);
    return gamma;
}

PV2D RK4Pusher2D::dudtAt(const PV2D & pos)
{
    const double uTheta = uThetaAt(pos);
    const double vTheta = uTheta / gammaAt(pos);
    const double z = pos.z;
    const double r = pos.r;
    const double ez = efield_->ez(z, r);
    const double er = efield_->er(z, r);
    return PV2D(
        -vTheta * magfield_->br(pos.z, pos.r) + Q0 / M0 * ez
        ,
        vTheta*(uTheta / r + Q0 / M0 * magfield_->bz(z, r)) + Q0 / M0 * er
    );
}

double RK4Pusher2D::uThetaAt(const PV2D & pos)
{
    return pTheta_/M0/pos.r - Q0 / M0 * magfield_->aTheta(pos.z, pos.r);
}
