#include <pusher/pusher.h>
#include <iostream>
#include <algorithm>

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
    auto aphi = APhiPusher(ef, mf);
    auto boris = BorisPusher(ef, mf);
    auto lf = LeapFrogPusher(ef, mf);
    auto rk = RK4Pusher(ef, mf);

#ifdef DEMO
    std::cout << "# 1:omega*dt 2:turns 3:rel_error_r_leapfrog 4:rel_error_r_boris 5:rel_error_r_aphi 6:rel_error_r_rk4\n";
#else
    std::cout << "# 1:omega*dt 2:tol_leapfrog 3:tol_boris 4:tol_aphi 5:tol_rk4\n";
#endif

    for (double dt = 1e-12; dt > 1e-18; dt *= 0.5) {
        const int64_t steps = static_cast<int64_t>(std::ceil(10*2*M_PI / (omega*dt)));
        std::clog << "omega*dt: " << omega*dt << '\n';
        const double startAngle = std::atan(dt * omega/2);
        aphi.setElectronInfo(0, r, 0, 0, aphi.pTheta(0, r, u), gamma);
        boris.setElectronInfo(r*std::cos(startAngle), r*std::sin(startAngle),0, 0, u, 0);
        lf.setElectronInfo(r*std::cos(startAngle), r*std::sin(startAngle),0, 0, u, 0);
        rk.setElectronInfo(r*std::cos(0), r*std::sin(0),0, 0, u, 0);
#ifdef DEMO
        if (dt < 1e-15)
            break;
        int64_t trigger = 1000;
        for (int64_t i = 1; i <= steps; ++i) {
            ++trigger;
            if (trigger*omega*dt > M_PI/60) { // every ~3 degree
                trigger = 0;
                const auto pa = aphi.pos();
                const auto pb = boris.pos();
                const auto pl = lf.pos();
                const auto pr = rk.pos();
                std::cout << omega*dt << ' ' << omega*dt*i/(2*M_PI) << ' '
                          <<  (fromPV3D(pl).r-r)/r << ' '
                          <<  (fromPV3D(pb).r-r)/r << ' '
                          <<  (pa.r-r)/r << ' '
                          <<  (fromPV3D(pr).r-r)/r << '\n';
            }
            aphi.step(dt);
            boris.step(dt);
            lf.step(dt);
            rk.step(dt);
        }
        std::cout << '\n';
#else
        double tola = 0;
        double tolb = 0;
        double toll = 0;
        double tolr = 0;
        for (int64_t i = 0; i < steps; ++i) {
            aphi.step(dt);
            boris.step(dt);
            lf.step(dt);
            rk.step(dt);
            const auto pa = aphi.pos();
            const auto pb = boris.pos();
            const auto pl = lf.pos();
            const auto pr = rk.pos();
            toll = std::max(toll, std::abs((fromPV3D(pl).r-r)/r));
            tolb = std::max(tolb, std::abs((fromPV3D(pb).r-r)/r));
            tola = std::max(tola, std::abs((pa.r-r)/r));
            tolr = std::max(tolr, std::abs((fromPV3D(pr).r-r)/r));
        }
        std::cout << omega*dt << ' '
                  << toll << ' '
                  << tolb << ' '
                  << tola << ' '
                  << tolr << '\n';
#endif
    }
    return 0;
}
