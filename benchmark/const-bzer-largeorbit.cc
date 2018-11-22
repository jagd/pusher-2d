#include <pusher/pusher.h>
#include <iostream>
#include <cmath>
#include <iomanip>


int main()
{
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    const double Bz = 1.0; // f ~ 28 GHz
    const double Er = -10e3 / 10e-2;
    const double gamma = 1.2;
    const double r = 0.001131326056780128; // Mathematica
    const double v = gamma2v(gamma);
    const double omega = v/r;
    const double u = v*gamma;

    std::clog << "r = " << r << '\n'
              << "v = " << v << '\n'
              << "gamma = " << v2gamma(v) << '\n'
              << "freq = " << omega/2/M_PI/1e9 << " GHz\n";

    const auto ef = std::make_shared<ConstErField>(Er);
    const auto mf = std::make_shared<ConstBzField>(Bz);
    auto boris = BorisPusher(ef, mf);
    auto lf = LeapFrogPusher(ef, mf);
    auto rk = RK4Pusher(ef, mf);
    auto aphi = APhiPusher(ef, mf);

    const double distance = std::abs((2*M_PI*10) / omega * v);
    const double totalEnergy = Q0*ef->pot(0, r) + (gamma-1)*(M0*C0*C0);

    for (double dt = 1e-12; dt > 1e-15; dt *= 0.5) {
        const int64_t steps = static_cast<int64_t>(distance / (v*dt));
        const double startAngle = std::atan(dt * omega / 2);
        boris.setElectronInfo(r*std::cos(startAngle), r*std::sin(startAngle), 0, 0, u, 0);
        lf.setElectronInfo(r*std::cos(startAngle), r*std::sin(startAngle), 0, 0, u, 0);
        rk.setElectronInfo(r*std::cos(0), r*std::sin(0), 0, 0, u, 0);
        aphi.setElectronInfo(0, r, 0, 0, aphi.pTheta(0, r, u), gamma);
        double trigger = 0;
#ifdef DEMO
        std::clog << "omega*dt = " << omega * dt << '\n';
#else
        std::cout << omega * dt << ' ';
#endif
        for (int i = 1; i <= steps; ++i) {
            aphi.step(dt);
            boris.step(dt);
            lf.step(dt);
            rk.step(dt);
            ++trigger;
#ifndef DEMO
        }
#else
            if (trigger*omega*dt > 0.1745) { // every ~10 degree
#endif
                const double potEnergyBoris = Q0 * ef->pot(boris.pos().z, fromPV3D(boris.pos()).r);
                const double kinEnergyBoris = (boris.gammaCurrent() - 1.0)*(M0*C0*C0);
                const double potEnergyRK = Q0 * ef->pot(rk.pos().z, fromPV3D(rk.pos()).r);
                const double kinEnergyRK = (rk.gammaCurrent() - 1.0)*(M0*C0*C0);
#ifdef DEMO
                std::cout << i * dt*omega / (2 * M_PI) << ' ';
#endif
                    std::cout
                    << fromPV3D(lf.pos()).r << ' '
                    << fromPV3D(boris.pos()).r << ' '
                    << aphi.pos().r << ' '
                    << fromPV3D(rk.pos()).r << ' '
                    << (potEnergyBoris + kinEnergyBoris - totalEnergy) / totalEnergy << ' '
                    << (potEnergyRK + kinEnergyRK - totalEnergy) / totalEnergy << '\n';
#ifdef DEMO
                trigger = 0;
            }
        }
        std::cout << '\n';
#endif
    }
    return 0;
}
