#include <pusher/pusher.h>
#include <iostream>
#include <algorithm>

int main()
{
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    const double B = 1; // f ~ 28 GHz
    const double u = 1.8e8; // 90 keV
    const double rLarmor = -M0/Q0/B * u;
    const double gamma = u2gamma(u);
    const double v = u/gamma;
    const double omega = v/rLarmor;

    const auto ef = std::make_shared<ConstEzField>(0);
    const auto mf = std::make_shared<ConstBzField>(B);
    auto lf2d = LeapFrogPusher2D(ef, mf);
    auto boris = BorisPusher(ef, mf);
    auto lf = LeapFrogPusher(ef, mf);
    auto rk = RK4Pusher(ef, mf);

    const double offset = 10000*rLarmor;

#ifdef DEMO
    std::cout << "# 1:omega*dt 2:turns 3:rel_error_r_leapfrog 4:rel_error_r_boris 5:rel_error_r_aphi 6:rel_error_r_rk4\n";
#else
    std::cout << "# 1:omega*dt 2:tol_r_leapfrog 3:tol_r_boris 4:tol_r_phi 5:tol_r_rk4 "
                 "6:tol_gamma_leapfrog  7:tol_gamma_boris 8:tol:gamma_aphi 9:tol_gamma_rk4\n";
#endif

    for (double dt = 1.3e-12; dt > 1e-14; dt *= std::pow(0.5, 1.0/8)) {
        const int64_t steps = static_cast<int64_t>(std::round(100*2*M_PI / (omega*dt)));
        std::clog << "omega*dt: " << omega*dt << '\n';
        const double startAngle = dt * omega/2;
        lf2d.setElectronInfo(0, std::sqrt(offset*offset +rLarmor*rLarmor + 2*offset*rLarmor*std::cos(startAngle)), 0, 0, lf2d.pTheta(0, offset+rLarmor, u), gamma);
        boris.setElectronInfo(offset+rLarmor*std::cos(startAngle), rLarmor*std::sin(startAngle),0, 0, u, 0);
        lf.setElectronInfo(offset+rLarmor*std::cos(startAngle), rLarmor*std::sin(startAngle),0, 0, u, 0);
        rk.setElectronInfo(offset+rLarmor*std::cos(0), rLarmor*std::sin(0),0, 0, u, 0);
        // in order to start from the same position
        for (int i = 0; i < 1024; ++i) {
            rk.step(dt / 2048);
        }
#ifdef DEMO
        if (dt < 1e-15)
            break;
        int64_t trigger = 1000;
        for (int64_t i = 1; i <= steps; ++i) {
            ++trigger;
            if (trigger*omega*dt > M_PI/60) { // every ~3 degree
                trigger = 0;
                const auto pl2 = lf2d.pos();
                const auto pb = boris.pos();
                const auto pl = lf.pos();
                const auto pr = rk.pos();
                // omega*dt*i, i start from 1 instead of 0,
                // because the first step was considered in the intialization.
                const double rSoll = std::sqrt(offset*offset +rLarmor*rLarmor + 2*offset*rLarmor*std::cos(omega*dt*i+startAngle));
                std::cout << omega*dt << ' '
                          << omega*dt*i/(2*M_PI) << ' '
                          << (fromPV3D(pl).r - rSoll) / rSoll << ' '
                          << (fromPV3D(pb).r-rSoll)/rSoll << ' '
                          << (pl2.r-rSoll)/rSoll << ' '
                          << (fromPV3D(pr).r - rSoll) / rSoll << '\n';
            }
            lf2d.step(dt);
            boris.step(dt);
            lf.step(dt);
            rk.step(dt);
        }
        std::cout << '\n';
#else
        double tolL2 = 0;
        double tolR2 = 0;
        double tolB = 0;
        double tolL = 0;
        double tolR = 0;
        double tolGL = 0;
        double tolGB = 0;
        double tolGA = 0;
        double tolGR = 0;
        for (int64_t i = 1; i <= steps; ++i) {
            lf2d.step(dt);
            boris.step(dt);
            lf.step(dt);
            rk.step(dt);
            // omega*dt*i, i start from 1 instead of 0,
            // because the first step was considered in the intialization.
            const double rSoll = std::sqrt(offset*offset +rLarmor*rLarmor + 2*offset*rLarmor*std::cos(omega*dt*i+startAngle));
            const auto pl2 = lf2d.pos();
            const auto pb = boris.pos();
            const auto pl = lf.pos();
            const auto pr = rk.pos();
            tolL2 = std::max(tolL2, std::abs((pl2.r-rSoll)/rSoll));
            tolB = std::max(tolB, std::abs((fromPV3D(pb).r-rSoll)/rSoll));
            tolL = std::max(tolL, std::abs((fromPV3D(pl).r-rSoll)/rSoll));
            tolR = std::max(tolR, std::abs((fromPV3D(pr).r-rSoll)/rSoll));
            tolGL = std::max(tolGL, std::abs(lf.gammaLastHalf() - gamma));
            tolGB = std::max(tolGB, std::abs(boris.gammaCurrent() - gamma));
            tolGA = std::max(tolGA, std::abs(lf2d.gammaCurrent() - gamma));
            tolGR = std::max(tolGR, std::abs(rk.gammaCurrent() - gamma));
        }
        std::cout << omega*dt << ' '
                  << tolL << ' '
                  << tolB << ' '
                  << tolL2 << ' '
                  << tolR << ' '
                  << tolGL << ' '
                  << tolGB << ' '
                  << tolGA << ' '
                  << tolGR << '\n';
#endif
    }
    return 0;
}
