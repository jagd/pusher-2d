#include <pusher/pusher.h>
#include <iostream>
#include <algorithm>

int main()
{
    const double B = 1.0; // f ~ 28 GHz
    const double u = 1e8; // f ~ 26 GHz considering gamma
    const double rLarmor = -M0/Q0/B * u;
    const double gamma = u2gamma(u);
    const double v = u/gamma;
    const double omega = v/rLarmor;

    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<ConstBzField>(B);
    auto aphi = APhiPusher(ef, mf);
    auto boris = BorisPusher(ef, mf);

    const double offset = 2*rLarmor;

#ifdef DEMO
    std::cout << "# 1:omega*dt 2:turns 3:error_r_boris 4:error_r_aphi\n";
#else
    std::cout << "# 1:omega*dt 2:tol_boris 3:tol_aphi\n";
#endif

    for (double dt = 1e-12; dt > 1e-16; dt *= 0.5) {
        const int64_t steps = static_cast<int64_t>(std::ceil(100*2*M_PI / (omega*dt)));
        std::clog << "omega*dt: " << omega*dt << '\n';
        const double startAngle = dt * omega/2;
        aphi.setElectronInfo(0, std::sqrt(offset*offset +rLarmor*rLarmor + 2*offset*rLarmor*std::cos(startAngle)), 0, 0, aphi.pTheta(0, offset+rLarmor, u), gamma);
        boris.setElectronInfo(offset+rLarmor*std::cos(startAngle), rLarmor*std::sin(startAngle),0, 0, u, 0);
#ifdef DEMO
        int64_t trigger = 1000;
        for (int64_t i = 1; i <= steps; ++i) {
            if (trigger*omega*dt > M_PI/60) { // every ~3 degree
                const auto pa = aphi.pos();
                const auto pb = boris.pos();
                // omega*dt*i, i start from 1 instead of 0,
                // because the first step was considered in the intialization.
                const double rSoll = std::sqrt(offset*offset +rLarmor*rLarmor + 2*offset*rLarmor*std::cos(omega*dt*i+startAngle));
                std::cout << omega*dt << ' '
                          << omega*dt*i/(2*M_PI) << ' '
                          <<  (fromPV3D(pb).r-rSoll)/rSoll << ' '
                          <<  (pa.r-rSoll)/rSoll << '\n';
                trigger = 0;
            }
            aphi.step(dt);
            boris.step(dt);
            ++trigger;
        }
        std::cout << '\n';
#else
        double tola = 0;
        double tolb = 0;
        for (int64_t i = 1; i <= steps; ++i) {
            aphi.step(dt);
            boris.step(dt);
            // omega*dt*i, i start from 1 instead of 0,
            // because the first step was considered in the intialization.
            const double rSoll = std::sqrt(offset*offset +rLarmor*rLarmor + 2*offset*rLarmor*std::cos(omega*dt*i+startAngle));
            const auto pa = aphi.pos();
            const auto pb = boris.pos();
            tola = std::max(tola, std::abs((pa.r-rSoll)/rSoll));
            tolb = std::max(tolb, std::abs((fromPV3D(pb).r-rSoll)/rSoll));
        }
        std::cout << omega*dt << ' ' << tolb << ' ' << tola << '\n';
#endif
    }
    return 0;
}
