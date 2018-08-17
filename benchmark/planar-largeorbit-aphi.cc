#include <pusher/pusher.h>
#include <iostream>

int main()
{
    const double B = 1.0; // f ~ 28 GHz
    const double u = 1e8; // f ~ 26 GHz considering gamma
    const double r = -M0/Q0/B * u;
    const double gamma = u2gamma(u);
    const double v = u/gamma;
    const double omega = v/r;

    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<ConstBzField>(B);
    auto pusher = APhiPusher(ef, mf);

    const double distance = v*1e-12 * 1000;

    for (double dt = 1e-11; dt > 1e-16; dt *= 0.8) {
        const int64_t steps = std::ceil(distance / (v*dt));
        std::clog << "steps: " << steps << " ; step width: " << dt*v << '\n';
        pusher.setElectronInfo(0, r, 0, 0, pusher.pTheta(0, r, u), gamma);
#ifdef DEMO
        int64_t trigger = 0;
        for (int64_t i = 0; i < steps; ++i) {
            pusher.step(dt);
            ++trigger;
            if (trigger*omega*dt > 0.1745) { // every ~10 degree
                const auto p = pusher.pos();
                std::cout << omega*dt*i << ' ' <<  (p.r-r)/r << '\n';
                trigger = 0;
            }
        }
        std::cout << '\n';
#else
        double tol = 0;
        for (int64_t i = 0; i < steps; ++i) {
            pusher.step(dt);
            const auto p = pusher.pos();
            tol = std::max(tol, std::abs((p.r-r)/r));
        }
        std::cout << omega*dt << ' ' << tol << '\n';
#endif
    }
    return 0;
}
