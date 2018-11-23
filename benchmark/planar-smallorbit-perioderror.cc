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

    const auto ef = std::make_shared<ConstEzField>(0);
    const auto mf = std::make_shared<ConstBzField>(B);
    auto aphi = APhiPusher(ef, mf);
    auto boris = BorisPusher(ef, mf);
    auto lf = LeapFrogPusher(ef, mf);
    auto rk = RK4Pusher(ef, mf);

    const double offset = 2*rLarmor;

    const double dt = 10/180.0*M_PI/omega;
    std::clog << "omega*dt = " << dt*omega <<  " ~ " << dt*omega*180/M_PI << " degree\n";
    const double startAngle = dt * omega/2;
    const int64_t steps = static_cast<int64_t>(std::round(30*2*M_PI/(omega*dt)));
    auto rAPhi = std::sqrt(offset*offset + rLarmor * rLarmor + 2 * offset*rLarmor*std::cos(startAngle));
    aphi.setElectronInfo(0, rAPhi, 0, 0, aphi.pTheta(0, offset+rLarmor, u), gamma);
    boris.setElectronInfo(offset+rLarmor*std::cos(startAngle), rLarmor*std::sin(startAngle),0, 0, u, 0);
    lf.setElectronInfo(offset+rLarmor*std::cos(startAngle), rLarmor*std::sin(startAngle),0, 0, u, 0);
    rk.setElectronInfo(offset+rLarmor*std::cos(0), rLarmor*std::sin(0),0, 0, u, 0);
    for (int i = 0; i < 1024; ++i) {
        rk.step(dt / 2048);
    }
    for (int64_t i = 1; i <= steps; ++i) {
        std::cout
            << omega * dt*i / (2 * M_PI) << ' '
            << fromPV3D(lf.pos()).r << ' '
            << fromPV3D(boris.pos()).r << ' '
            << aphi.pos().r << ' '
            << fromPV3D(rk.pos()).r << '\n';
        aphi.step(dt);
        boris.step(dt);
        lf.step(dt);
        rk.step(dt);
    }
    return 0;
}
