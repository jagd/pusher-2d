#include <pusher/pusher.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>


int main()
{
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    const double r = 1;
    const double Ez = -1e8;
    const double uz0 = -1e8; // free to change
    const double exactETotal = (std::sqrt((uz0*uz0 /C0/C0)+1)-1) *M0*C0*C0 / Q0 * 1e-3;
    const auto ef = std::make_shared<ConstEzField>(Ez);
    const auto mf = std::make_shared<ConstBzField>(0);

    auto lf = LeapFrogPusher(ef, mf);
    auto boris = BorisPusher(ef, mf);
    auto aphi = APhiPusher(ef, mf);
    auto rk = RK4Pusher(ef, mf);

#ifdef DEMO
    const double dtScale = 256;
#else
    for (double dtScale = 8192*2; dtScale > 0.01; dtScale /= 2) {
#endif
        const int64_t steps = static_cast<int64_t>((1 << 14) / dtScale);
        const double dt = 1e-15*dtScale;
        const double uInitHalf = uz0 - Q0 * Ez / M0 * dt / 2;
        lf.setElectronInfo(r, 0, 0, 0, 0, uInitHalf);
        boris.setElectronInfo(r, 0, 0, 0, 0, uInitHalf);
        aphi.setElectronInfo(0, r, uInitHalf, 0, aphi.pTheta(0, r, 0), u2gamma(uz0));
        rk.setElectronInfo(r, 0, 0, 0, 0, uz0);
        std::clog << "dt = " << dt << '\n';
        int64_t j = 0;
        for (int64_t i = 1; i <= steps; ++i) {
            lf.step(dt);
            boris.step(dt);
            aphi.step(dt);
            rk.step(dt);
#ifndef DEMO
        }
        const int64_t i = steps;
        j = static_cast<int64_t>(128 / dtScale);
#endif
        if (++j >= (128 / dtScale)) {
            j = 0;
            std::clog << i * 100.0 / steps << " %\n" << aphi.gammaCurrent() << '\n';
            const double potEnergyBoris = Q0 * ef->pot(boris.pos().z, 0);
            const double kinEnergyBoris = (boris.gammaCurrent() - 1.0)*(M0*C0*C0);
            const double potEnergyRK = Q0 * ef->pot(rk.pos().z, 0);
            const double kinEnergyRK = (rk.gammaCurrent() - 1.0)*(M0*C0*C0);
            const double t = dt * (i);
            const double exactZ = C0 /Q0 /Ez*(std::sqrt(M0*M0*C0*C0 + std::pow(Q0*Ez*t+uz0*M0, 2)) -std::sqrt(M0*M0*C0*C0+uz0*uz0*M0*M0));
            std::cout
                << dt << ' '
                << i * dt << ' '
                << exactZ << ' '
                << lf.pos().z << ' '
                << boris.pos().z << ' '
                << aphi.pos().z << ' '
                << rk.pos().z << ' '
                << exactETotal << ' '
                << (potEnergyBoris + kinEnergyBoris) / Q0 / 1e3 << ' '
                << (potEnergyRK + kinEnergyRK) / Q0 / 1e3 << '\n';
        }
#ifdef DEMO
        }
#endif
#ifndef DEMO
    } // for dt
#endif
    return 0;
}
