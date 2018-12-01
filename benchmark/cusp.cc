#include <pusher/pusher.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>

int main()
{
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    const double b1 = -0.01;
    const double b2 = -b1;
    const double z0field = 1e-4;
    const double r0 = 1e-2;
    const double u0 = 1e8;
    const double z0 = -0.1;

    const auto ef = std::make_shared<ConstEzField>(0);
    const auto mf = std::make_shared<BiUniformMagField>(z0field, b1, b2);
    auto aphi = LeapFrogPusher2D(ef, mf);
    auto boris = BorisPusher(ef, mf);
    auto rk = RK4Pusher(ef, mf);
    auto lf = LeapFrogPusher(ef, mf);

    // const double flux0 = r0 * r0*mf->bz(0, r0);

#ifdef DEMO
    const double dtScale = 1;
#else
    for (double dtScale = 2; dtScale > 1e-3; dtScale /= 2) {
        double errL = 0;
        double errB = 0;
        double errA = 0;
        double errR = 0;
        double stepL = 0;
        double stepB = 0;
        double stepA = 0;
        double stepR = 0;
#endif
        const double dt = 1e-12 * dtScale;
        const int64_t maxStep = static_cast<int64_t>((1 << 13) / dtScale);
        std::clog << "dt = " << dt << '\n';

        lf.setElectronInfo(r0, 0, z0, 0, 0, u0);
        boris.setElectronInfo(r0, 0, z0, 0, 0, u0);
        rk.setElectronInfo(r0, 0, z0, 0, 0, u0);
        aphi.setElectronInfo(z0, r0, u0, 0, aphi.pTheta(z0, r0, 0), u2gamma(u0));

        for (int64_t i = 0; i < maxStep; ++i) {
#ifdef DEMO
            if (i % static_cast<int64_t>(std::ceil(64 / dtScale)) == 0) {
                std::clog << i * 100.0 / maxStep << " %\n" << '\n';
                std::cout
                    << i * dt << ' '
                    << lf.pos().z << ' '
                    << fromPV3D(lf.pos()).r << ' '
                    << boris.pos().z << ' '
                    << fromPV3D(boris.pos()).r << ' '
                    << aphi.pos().z << ' '
                    << aphi.pos().r << ' '
                    << rk.pos().z << ' '
                    << fromPV3D(rk.pos()).r << '\n';
            }
#else
            errL = std::max(errL, std::abs(fromPV3D(lf.pos()).r - r0));
            errB = std::max(errB, std::abs(fromPV3D(boris.pos()).r - r0));
            errA = std::max(errA, std::abs(aphi.pos().r - r0));
            errR = std::max(errR, std::abs(fromPV3D(rk.pos()).r - r0));
            if (std::abs(lf.pos().z) < z0field)
                ++stepL;
            if (std::abs(boris.pos().z) < z0field)
                ++stepB;
            if (std::abs(aphi.pos().z) < z0field)
                ++stepA;
            if (std::abs(rk.pos().z) < z0field)
                ++stepR;
#endif
            aphi.step(dt);
            boris.step(dt);
            rk.step(dt);
            lf.step(dt);
        }
#ifndef DEMO
        std::cout
            << dt << ' '
            << stepL << ' ' << errL/r0 << ' '
            << stepB << ' ' << errB/r0 << ' '
            << stepA << ' ' << errA/r0 << ' '
            << stepR << ' ' << errR/r0 << '\n';
    }
#endif
    return 0;
}
