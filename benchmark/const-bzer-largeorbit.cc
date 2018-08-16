#include <pusher/pusher.h>
#include <iostream>
#include <cmath>
#include <iomanip>


int main()
{
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
    auto aphi = APhiPusher(ef, mf);

    const double distance = std::abs((2*M_PI*10) / omega * v);
    const double totalEnergy = Q0*ef->pot(0, r) + (gamma-1)*(M0*C0*C0);

    for (double dt= 1e-12; dt > 1e-15; dt *= 0.8) {
        const int64_t steps = static_cast<int64_t>(std::ceil(distance / (v*dt)));
        const double startAngle = std::atan(dt * omega/2);
        boris.setElectronInfo(r*std::cos(startAngle), r*std::sin(startAngle),0, 0, u, 0);
        aphi.setElectronInfo(0, r, 0, 0, aphi.pTheta(0, r, u));
        double trigger = 0;
        std::clog << "dt = " << dt << '\n';
        for (int i = 1; i <= steps; ++i) {
            aphi.step(dt);
            boris.step(dt);
            ++trigger;
            if (trigger*omega*dt > 0.1745) { // every ~10 degree
                const auto pBoris = boris.pos();
                const auto pAPhi = aphi.pos();
                const double potEnergyBoris = Q0*ef->pot(pBoris.z, 0);
                const double kinEnergyBoris = (boris.gamma()-1.0)*(M0*C0*C0);
                std::cout << std::setprecision(18)
                          // Col 1: time
                          << i*dt << ' '
                          // Col 2-3: Radius of both methods
                          <<  std::sqrt(pBoris.x*pBoris.x + pBoris.y*pBoris.y) << ' '
                          <<  pAPhi.r << ' '
                          // Col 4-5: Relative error of radius
                          << (std::sqrt(pBoris.x*pBoris.x + pBoris.y*pBoris.y) - r)/r << ' '
                          << (pAPhi.r - r)/r  << ' '
                          // Col 6: Error of energy
                          << potEnergyBoris + kinEnergyBoris - totalEnergy / totalEnergy << '\n';
                trigger = 0;
            }
        }
        std::cout << '\n';
    }
    return 0;
}
