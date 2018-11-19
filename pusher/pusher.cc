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
    //
    // Example for an electron (Q > 0):
    // phi increases -> Epot=phi*Q decreases -> Ekin=Etotal-Epot increases
    //               -> Gamma=Ekin/(M0C^2) increases -> dg/dphi positive
    // For a positive ion (Q < 0):
    // phi increases -> Epot=phi*Q increases -> Ekin=Etotal-Epot decreases
    //               -> Gamma=Ekin/(M0C^2) decreases -> dg/dphi negative
    // Hence, the sign of dgamma/dphi is in -1*Q0.
    //
    const double z = pos_.z;
    const double r = pos_.r;
    const double p = pTheta_;
    const static double dgdphi = -Q0/(M0*C0*C0);
    const double uTheta = pTheta_ / (M0*r) - Q0 / M0 * magfield_->aTheta(z, r);
    const double commonTerm = g*(C0*C0);
    const double dgdz = -dgdphi*efield_->ez(z, r);
    const double dgdr = -dgdphi*efield_->er(z, r);
    const PV2D gradient(
    	Q0 / M0 * uTheta*magfield_->aTheta2z(z, r) + dgdz*commonTerm
   	,
    	Q0 / M0 * uTheta*magfield_->aTheta2r(z, r) + pTheta_ * uTheta / (r*r*M0) + dgdr*commonTerm
    );
    const PV2D uNextHalf = uLastHalf_ + dt * gradient / g;
#ifdef GAMMA_CORRECTION_APHI
    const double disc = g * g + 2 * dt* (dgdz*uNextHalf.z + dgdr*uNextHalf.r);
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
	uLast_ = PV3D(ux, uy, uz);
}

static PV3D p2v(const PV3D &p)
{
	const auto u = p / M0 / C0;
	const auto gamma = std::sqrt(1 + norm2(u));
	return p / M0 / gamma;
}

void RK4Pusher::step(double dt)
{
	const double gammaLast = usqr2gamma(norm2(uLast_));
	const PV3D vLast = uLast_ / gammaLast;
	const PV3D pLast = uLast_ * M0;
	
	const PV3D r0 = dt * p2v(pLast);
	const PV3D p0 = dt * f3D(pos_, pLast);

	const auto r1 = dt * p2v(pLast + 0.5*p0);
	const auto p1 = dt * f3D(pos_ + 0.5*r0, pLast + 0.5*p0);
	const auto r2 = dt * p2v(pLast + 0.5*p1);
	const auto p2 = dt * f3D(pos_ + 0.5*r1, pLast + 0.5*p1);
	const auto r3 = dt * p2v(pLast + p2);
	const auto p3 = dt * f3D(pos_ + r2, pLast + p2);

	const PV3D pNext = pLast + 1.0 / 6 * (p0 + 2 * p1 + 2 * p2 + p3);
	const PV3D vNext = p2v(pNext);
	pos_ += 1.0 / 6 * (r0 + 2 * r1 + 2 * r2 + r3);
	uLast_ = vNext * v2gamma(std::sqrt(norm2(vNext)));
}

PV3D RK4Pusher::pos() const
{
	return pos_;
}

PV3D RK4Pusher::uCurrent() const
{
	return uLast_;
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

PV3D RK4Pusher::f3D(const PV3D & pos, const PV3D & p) const
{
	return Q0 * (cross(p2v(p), b3D(pos)) + e3D(pos));
}

double RK4Pusher::gammaCurrent() const
{
    return usqr2gamma(norm2(uLast_));
}
