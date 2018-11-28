#include <pusher/pusher.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>


int main()
{
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    const double r = 1;
    const double Ez0 = -1e8;
    const double Ek = 1e6;
    const double uz0 = -1e8; // free to change
    const auto ef = std::make_shared<LinearEzField>(Ez0, Ek);
    const auto mf = std::make_shared<ConstBzField>(0);

    auto lf = LeapFrogPusher(ef, mf);
    auto boris = BorisPusher(ef, mf);
    auto aphi = APhiPusher(ef, mf);
    auto rk = RK4Pusher(ef, mf);

    for (double dtScale = 8192 * 2; dtScale > 0.01; dtScale /= 2) {
        const int64_t steps = static_cast<int64_t>((1 << 14) / dtScale);
        const double dt = 1e-15*dtScale;
        const double uInitHalf = uz0 - Q0 * Ez0 / M0 * dt / 2;
        lf.setElectronInfo(r, 0, 0, 0, 0, uInitHalf);
        boris.setElectronInfo(r, 0, 0, 0, 0, uInitHalf);
        aphi.setElectronInfo(0, r, uInitHalf, 0, aphi.pTheta(0, r, 0), u2gamma(uz0));
        rk.setElectronInfo(r, 0, 0, 0, 0, uz0);
        for (int64_t i = 1; i <= steps; ++i) {
            lf.step(dt);
            boris.step(dt);
            aphi.step(dt);
            rk.step(dt);
        }

        std::clog << "dt = " << dt << " gamma = " << aphi.gammaCurrent() << '\n';
        std::cout
            << dt << ' '
            << lf.pos().z << ' '
            << boris.pos().z << ' '
            << aphi.pos().z << ' '
            << rk.pos().z << '\n';
    }
    return 0;
}
