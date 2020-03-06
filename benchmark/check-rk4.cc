#include <pusher/pusher.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <limits>

#define C511 510.99895000

PV3D calcU(const double Ekin_energy, const double alpha_number)
{
    //calculate helper variables
    const auto beta0 = std::sqrt( 1.0 - (C511*C511) / std::pow(C511+Ekin_energy,2) );
    const auto gamma0_number = 1.0/std::sqrt(1.0-(beta0*beta0));
    const auto u0z_velocity = gamma0_number * C0 * std::sqrt( (beta0*beta0) / (alpha_number*alpha_number + 1.0) );
    const auto u0t_velocity = gamma0_number * C0 * std::sqrt( (beta0*beta0) / (1.0/(alpha_number*alpha_number) + 1.0) );

    const PV3D u(u0t_velocity, 0, u0z_velocity);
    return u;
}

double lamourRadius(const double _B_bfield, const double _Ekin_energy, const double _alpha)
{
    const double beta0 = std::sqrt( 1.0 - (C511*C511) / std::pow(C511+_Ekin_energy,2) );
    const double gamma0_number = 1.0/std::sqrt(1.0-(beta0*beta0));
    const double v0t_velocity = C0 * std::sqrt( (beta0*beta0) / (1.0/(_alpha*_alpha) + 1.0) );
    const double omegaH = (_B_bfield * std::abs(Q0)) / (gamma0_number * M0);

    return v0t_velocity / omegaH;
}


int main()
{
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    const double Bz = 0.9;
    const double dt = 0.1e-11;
    const double eKin = 30;
    const double alpha = 1.1;

    const auto ef = std::make_shared<ConstErField>(0);
    const auto mf = std::make_shared<ConstBzField>(Bz);

    auto rk = RK4Pusher(ef, mf);

    const auto u = calcU(eKin, alpha);
    const auto rLarmor = lamourRadius(Bz, eKin, alpha);

    rk.setElectronInfo(rLarmor, 0, 0, u.x, u.y, u.z);

    for (int i = 0; i <= 0.5e-9/dt; ++i) {
        const double kinEnergyRK = (rk.gammaCurrent() - 1.0)*(M0*C0*C0);
        std::cout
                << fromPV3D(rk.pos()).r
                << ' '
                << kinEnergyRK * 6.242e+18 / 1000
                << std::endl; // will be flushe
        rk.step(dt);
    }
    return 0;
}
