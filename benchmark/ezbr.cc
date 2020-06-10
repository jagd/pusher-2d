/**
  This benchmark considers a uniform axial electric field with grad(E) to be a constant (or zero).
  There can be a magnetic field with a constant Br, which can be zero, too.
*/

#include <pusher/pusher.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>


int main()
{
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    const double r = 1e9;
    const double Ez0 = -1e8;
    const double Ek = 0; // free to change, in the test case 1e6
    const double uz0 = -1e8; // free to change
    const double exactETotal = (std::sqrt((uz0*uz0 /C0/C0)+1)-1) *M0*C0*C0 / Q0 * 1e-3;
    const auto ef = std::make_shared<LinearEzField>(Ez0, Ek);
    // free to change, Br=0 will yield the exact solution, which is independent of initial velocity
    const auto mf = std::make_shared<ConstBrField>(0);

    auto lf = LeapFrogPusher(ef, mf);
    auto boris = BorisPusher(ef, mf);
    auto aphi = LeapFrogPusher2D(ef, mf);
    auto rk = RK4Pusher(ef, mf);

#ifdef DEMO
    const double dtScale = 256;
#else
    for (double dtScale = 8192*2; dtScale > 0.001; dtScale /= 2) {
#endif
        const int64_t steps = static_cast<int64_t>((1 << 14) / dtScale);
        const double dt = 1e-15*dtScale;
        PV3D uInitHalf(0, 0, uz0 - Q0 * Ez0 / M0 * dt / 2); // exact for Br = 0 && Ek == 0
        if (mf->br(0,r) != 0 || Ek != 0) {
#if 0
            const int N = 1000000;
#else
            // I fix the scenario-dependent number dt~8e-13 for better accuracy
            const int N = std::max((int)std::ceil(dt/8e-13), 2);
#endif
            const double reverseSubStep = -dt / 2 / N;
            rk.setElectronInfo(r, 0, 0, 0, 0, uz0);
            for (int i = 0; i < N; ++i) {
                rk.step(reverseSubStep);
            }
            uInitHalf = rk.uCurrent();
        }
        lf.setElectronInfo(r, 0, 0, uInitHalf.x, uInitHalf.y, uInitHalf.z);
        boris.setElectronInfo(r, 0, 0, uInitHalf.x, uInitHalf.y, uInitHalf.z);
        aphi.setElectronInfo(0, r, uInitHalf.z, -fromPV3D(uInitHalf).r, aphi.pTheta(0, r, 0), u2gamma(uz0));
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
            const double exactZ = C0 /Q0 /Ez0*(std::sqrt(M0*M0*C0*C0 + std::pow(Q0*Ez0*t+uz0*M0, 2)) -std::sqrt(M0*M0*C0*C0+uz0*uz0*M0*M0));
            std::cout
                << dt << ' '
                << i * dt << ' '
                << exactZ << ' ' // not exact anymore if Br != 0 !!!
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
