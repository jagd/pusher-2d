#include "pusher.h"
#include "pv.h"
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
}


void BorisPusher::step(double dt)
{
    const PV2D zr = fromPV3D(pos_);
    const PV3D posNorm = (zr.r == 0) ? PV3D(0,0,1) : (pos_/std::sqrt(dot(pos_, pos_)));
    const auto er = efield_->er(zr.z, zr.r);
    const PV3D e3d(er*posNorm.x, er*posNorm.y, efield_->ez(zr.z, zr.r));
    const PV3D uMinus = uLastHalf_ + (Q0/2/M0)*dt*e3d;
    const double uMinusNorm2 = dot(uMinus, uMinus);
    const double gammaN = std::sqrt(1.0 + uMinusNorm2/(C0*C0));

    const auto br = magfield_->br(zr.z, zr.r);
    const PV3D b3d(br*posNorm.x, br*posNorm.y, magfield_->bz(zr.z, zr.r));
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
    const double prqamg = (
        pTheta_/pos_.r - Q0*magfield_->aTheta(pos_.z, pos_.r)
    ) / (M0*g);
    //
    // Example for an electron (Q > 0):
    // phi increases -> Epot=phi*Q decreases -> Ekin=Etotal-Epot increases
    //               -> Gamma=Ekin/(M0C^2) increases -> dg/dphi positive
    // For a positive ion (Q < 0):
    // phi increases -> Epot=phi*Q increases -> Ekin=Etotal-Epot decreases
    //               -> Gamma=Ekin/(M0C^2) decreases -> dg/dphi negative
    // Hence, the sign of dgamma/dphi is in -1*Q0.
    //
    const double dgdphi = -Q0/(M0*C0*C0);
    const double commonTerm = (prqamg*prqamg*(M0/2)*dgdphi - Q0);
    const double gradFz = -Q0*prqamg*magfield_->aTheta2z(pos_.z, pos_.r)
        + commonTerm * efield_->ez(pos_.z, pos_.r);
    const double gradFr = (pos_.r == 0) ? 0 :
        -prqamg*(pTheta_/(pos_.r*pos_.r)+Q0*magfield_->aTheta2r(pos_.z, pos_.r))
        + commonTerm * efield_->er(pos_.z, pos_.r);
    const PV2D f(-gradFz, -gradFr);
    const PV2D uNextHalf = f/M0*dt + uLastHalf_;
#ifdef APHI_GAMMA_CORRECTION
    const PV2D halfDist = (uNextHalf+uLastHalf_)*(dt/4/g);
    // gamma_corrected_{t+dt/2} = gaemma_t
    //     + grad{gamma}|_{t=t} (dot) (u_{t-dt/2}+ u_{t+dt/2})/2/gamma_{t} * dt/2
    // The next minus is for grad(phi) = -E
    const double gammaNextHalf = g - (
        dgdphi*efield_->ez(pos_.z, pos_.r)*halfDist.z
      + dgdphi*efield_->er(pos_.z, pos_.r)*halfDist.r
    );
    pos_ = uNextHalf/gammaNextHalf*dt + pos_;
#else
    pos_ = uNextHalf/g*dt + pos_;
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
    pos_ = PV2D(z, r);
    uLastHalf_ = PV2D(uzLastHalf, urLastHalf);
    pTheta_ = pTheta;
    if (gamma < 1.0) {
        const double uTheta = (pTheta / r-(magfield_->aTheta(z, r)*Q0))/M0;
        gamma = u2gamma(std::sqrt(uzLastHalf*uzLastHalf + urLastHalf*urLastHalf + uTheta*uTheta));
    }
    totalEnergy_ = (gamma - 1) * (M0*C0*C0) + Q0*efield_->pot(z, r);
}

double APhiPusher::gamma() const
{
    const double Ekin = totalEnergy_-Q0*efield_->pot(pos_.z, pos_.r);
//    assert(Ekin >= 0);
    const double gamma = 1.0 + Ekin/(M0*C0*C0);
    return gamma;
}

PV2D APhiPusher::u() const
{
    return uLastHalf_;
}

PV2D APhiPusher::pos() const
{
    return pos_;
}
