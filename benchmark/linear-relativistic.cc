#include <pusher/pusher.h>
#include <iostream>
#include <iomanip>


int main()
{
    const double r = 1;
    const double Ez = -1e8;
    const auto ef = std::make_shared<ConstEzField>(Ez);
    const auto mf = std::make_shared<ConstBzField>(0);

    auto lf = LeapFrogPusher(ef, mf);
    auto boris = BorisPusher(ef, mf);
    auto aphi = APhiPusher(ef, mf);
    auto rk = RK4Pusher(ef, mf);

    const double dtScale = 4096;
    const int64_t steps = static_cast<int64_t>((1 << 14)/dtScale);
    const double dt = 1e-15*dtScale;
    const double vInitHalf = -Ez * Q0/M0*dt / 2 / (std::sqrt(1+(dt*dt/4)*Ez*Ez*Q0*Q0/(M0*M0*C0*C0)));
    const double uInitHalf = v2u(vInitHalf);
    lf.setElectronInfo(r, 0, 0, 0, 0, uInitHalf);
    boris.setElectronInfo(r, 0, 0, 0, 0, uInitHalf);
    aphi.setElectronInfo(0, r, uInitHalf, 0, aphi.pTheta(0, r, 0), u2gamma(0));
    rk.setElectronInfo(r, 0, 0, 0, 0, 0);
    std::clog << "dt = " << dt << '\n';
    int64_t j = 0;
    for (int64_t i = 1; i <= steps; ++i) {
        lf.step(dt);
        boris.step(dt);
        aphi.step(dt);
        rk.step(dt);
        if (++j >= (128/dtScale)) {
            j = 0;
            std::clog << i * 100.0 / steps << " %\n" << aphi.gammaCurrent() << '\n';
            const double potEnergyBoris = Q0 * ef->pot(boris.pos().z, 0);
            const double kinEnergyBoris = (boris.gammaCurrent() - 1.0)*(M0*C0*C0);
            const double potEnergyRK = Q0 * ef->pot(rk.pos().z, 0);
            const double kinEnergyRK = (rk.gammaCurrent() - 1.0)*(M0*C0*C0);
            const double t = dt * (i);
            const double exactZ = (M0*C0*C0/Ez/Q0)*(std::sqrt(1+(Ez*Ez*Q0*Q0/(M0*M0*C0*C0)*t*t)) - 1);
            std::cout << i * dt << ' '
                << exactZ << ' '
                << lf.pos().z << ' '
                << boris.pos().z << ' '
                << aphi.pos().z << ' '
                << rk.pos().z << ' '
                << (potEnergyBoris + kinEnergyBoris)/Q0/1e3<< ' '
                << (potEnergyRK + kinEnergyRK)/Q0/1e3<< '\n';
        }
    }
    return 0;
}
