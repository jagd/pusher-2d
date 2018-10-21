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
    const double uMinusNorm2 = dot(uMinus, uMinus);
    const double gammaN = std::sqrt(1.0 + uMinusNorm2/(C0*C0));

    const auto br = magfield_->br(zr.z, zr.r);
    const PV3D b3d(br*planarNorm.x, br*planarNorm.y, magfield_->bz(zr.z, zr.r));
    const PV3D t = (Q0/M0/2)*dt/gammaN*b3d;
    const PV3D s = 2.0/(1+dot(t,t)) * t;
    const PV3D uPrime = uMinus + cross(uMinus, t);
    const PV3D uPlus = uMinus + cross(uPrime, s);

    const PV3D uNextHalf = uPlus + (Q0/2/M0)*dt*e3d;
    pos_ = pos_ + dt/(std::sqrt(1+dot(uNextHalf, uNextHalf)/(C0*C0))) * uNextHalf;
    uLastHalf_ = uNextHalf;
    gammaAtPos_ = gammaN;
}


PV3D BorisPusher::pos() const
{
    return pos_;
}


PV3D BorisPusher::u() const
{
    return uLastHalf_;
}


double BorisPusher::gamma() const
{
    return gammaAtPos_;
}


double APhiPusher::pTheta(double z, double r, double uTheta) const
{
    return r*(uTheta*M0 + magfield_->aTheta(z, r)*Q0);
}


void APhiPusher::step(double dt)
{
    const double g = gamma();
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
    const double dgdphi = -Q0/(M0*C0*C0);
    const double mur = p - Q0*r*magfield_->aTheta(z, r); // m * u_theta * r
    const double commonTerm = p*p + C0*C0*M0*M0*r*r - Q0*r*magfield_->aTheta(z, r)*(p+mur);
    const double denom = M0*M0*r*r*g*g*g;
    // The minus in the middle is due to Grad(phi) = -E
    const double fz = (
        Q0*r*g*mur*magfield_->aTheta2z(z, r)
        - commonTerm*dgdphi*efield_->ez(z, r)
    ) / denom;
    const double fr = (
        g*mur*(p+Q0*r*r*magfield_->aTheta2r(z, r))
        - r*commonTerm*dgdphi*efield_->er(z, r)
    ) / (denom*r);
    const PV2D f(fz, fr);
    vLastHalf_ += f*dt;
    pos_ += vLastHalf_*dt;
}


void APhiPusher::setElectronInfo(
    double z,
    double r,
    double vzLastHalf,
    double vrLastHalf,
    double pTheta,
    double gamma
)
{
    assert(gamma >= 1);
    pos_ = PV2D(z, r);
    vLastHalf_ = PV2D(vzLastHalf, vrLastHalf);
    pTheta_ = pTheta;
    totalEnergy_ = (gamma - 1) * (M0*C0*C0) + Q0*efield_->pot(z, r);
}

double APhiPusher::gamma() const
{
    const double Ekin = totalEnergy_-Q0*efield_->pot(pos_.z, pos_.r);
//    assert(Ekin >= 0);
    const double gamma = 1.0 + Ekin/(M0*C0*C0);
    return gamma;
}

PV2D APhiPusher::v() const
{
    return vLastHalf_;
}

PV2D APhiPusher::pos() const
{
    return pos_;
}

void LeapFrog::setElectronInfo(double x, double y, double z, double ux, double uy, double uz)
{
	pos_ = PV3D(x, y, z);
	uLastHalf_ = PV3D(ux, uy, uz);
}

void LeapFrog::step(double dt)
{
	// not optimized, just for reference
    // hence, v and u are (unnecessary) converted back and forth 
    const PV2D zr = fromPV3D(pos_);
    const PV3D planarNorm = (zr.r == 0) ? PV3D(0,0,1) : (pos_/zr.r);
    const auto br = magfield_->br(zr.z, zr.r);
    const PV3D b3d(br*planarNorm.x, br*planarNorm.y, magfield_->bz(zr.z, zr.r));
    const auto er = efield_->er(zr.z, zr.r);
    const PV3D e3d(er*planarNorm.x, er*planarNorm.y, efield_->ez(zr.z, zr.r));
	const double gammaLastHalf = u2gamma(std::sqrt(dot(uLastHalf_, uLastHalf_)));
	const PV3D vLastHalf = uLastHalf_ / gammaLastHalf;
	const PV3D a = Q0/M0 * (cross(vLastHalf, b3d) + e3d);
	const PV3D vNextHalf_ = vLastHalf + a * dt;
	uLastHalf_ = vNextHalf_ * v2gamma(std::sqrt(dot(vNextHalf_, vNextHalf_)));
	pos_ += vNextHalf_ * dt;
}


PV3D LeapFrog::pos() const
{
	return pos_;
}

PV3D LeapFrog::u() const
{
	return uLastHalf_;
}
