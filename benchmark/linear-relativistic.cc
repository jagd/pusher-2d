#include <pusher/pusher.h>
#include <iostream>
#include <iomanip>


int main()
{
    const double zInit = 0;
    const double r = 1;
    const double uz0 = 0;
    const double ur0 = 1e7;
    const double u0 = std::sqrt(uz0*uz0 + ur0 * ur0);
    const auto ef = std::make_shared<ConstEzField>(-1000e3/10e-2);
    const auto mf = std::make_shared<ConstBzField>(0);
    auto boris = BorisPusher(ef, mf);

    auto aphi = APhiPusher(ef, mf);

    const double dtScale = 0.01;
    const int64_t steps = static_cast<int64_t>((1 << 22)/dtScale);
    double dt = 1e-15;
    boris.setElectronInfo(r, 0, zInit, ur0, 0, uz0);
    aphi.setElectronInfo(zInit, r, uz0, ur0, aphi.pTheta(zInit, r, 0), u2gamma(u0));
    std::clog << "dt = " << dt << '\n';
    for (int i = 0; i < steps; ++i) {
        aphi.step(dt);
        boris.step(dt);
        if (i % static_cast<int64_t>(std::ceil(256 / dtScale)) == 0) {
            std::clog << i * 100.0 / steps << " %\n";
            const auto pBoris = boris.pos();
            const auto pAPhi = aphi.pos();
            const double potEnergyBoris = Q0 * ef->pot(pBoris.z, 0);
            const double kinEnergyBoris = (boris.gammaCurrent() - 1.0)*(M0*C0*C0);
            std::cout << i * dt << ' '
                << boris.pos().x << ' '
                << boris.pos().z << ' '
                << aphi.pos().r << ' '
                << aphi.pos().z << ' '
                << potEnergyBoris + kinEnergyBoris<< '\n';
        }
    }
    return 0;
}
