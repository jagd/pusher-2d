#include <pusher/pusher.h>
#include <iostream>

int main()
{
    const double B = 1.0; // f ~ 28 GHz
    const double u = 1e8; // f ~ 26 GHz considering gamma
    const double rLarmor = -M0/Q0/B * u;
    const double gamma = u2gamma(u);
    const double v = u/gamma;
    const double omega = v/rLarmor;

    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<HomogeneousMagField>(B);
    auto pusher = BorisPusher(ef, mf);

    const double distance = v*1e-12 * 1000;
    const double offset = 2*rLarmor;

    for (double dt = 1e-11; dt > 1e-16; dt *= 0.8) {
        const int64_t steps = static_cast<int64_t>(std::ceil(distance / (v*dt)));
        std::clog << "steps: " << steps << " ; 3D step width: " << dt*v << '\n';
        const double startAngle = dt * omega/2;
        pusher.setElectronInfo(offset+rLarmor*std::cos(startAngle), rLarmor*std::sin(startAngle),0, 0, u, 0);
#ifdef DEMO
        int64_t trigger = static_cast<int64_t>(6.0/omega/dt);
        for (int64_t i = 0; i < steps; ++i) {
            if (trigger*omega*dt > 0.1745) { // every ~10 degree
                const auto p = pusher.pos();
                std::cout << omega*dt*i << ' ' << std::sqrt(p.x*p.x+p.y*p.y) << '\n';
                trigger = 0;
            }
            pusher.step(dt);
            ++trigger;
        }
        std::cout << '\n';
#else
        double tol = 0;
        for (int64_t i = 0; i < steps; ++i) {
            pusher.step(dt);
            const auto p = pusher.pos();
            const double rSoll = std::sqrt(offset*offset +rLarmor*rLarmor + 2*offset*rLarmor*std::cos(startAngle+omega*dt*(i+1)));
            tol = std::max(tol, std::abs((std::sqrt(p.x*p.x+p.y*p.y)-rSoll)/rSoll));
        }
        std::cout << omega*dt << ' ' << tol << '\n';
#endif
    }
    return 0;
}
