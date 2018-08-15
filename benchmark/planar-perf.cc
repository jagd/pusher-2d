#include <pusher/pusher.h>
#include <iostream>
#include <chrono>

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

    const double offset = 2*rLarmor;
    const int64_t steps = 20000000;
    const double dt = 1e-15;
    const double startAngle = dt * omega/2;

#ifndef BORIS_ONLY
    auto aphi = APhiPusher(ef, mf);
    aphi.setElectronInfo(0, std::sqrt(offset*offset +rLarmor*rLarmor + 2*offset*rLarmor*std::cos(startAngle)), 0, 0, aphi.pTheta(0, offset+rLarmor, u));
    const auto aphiBegin = std::chrono::steady_clock::now();
    for (int64_t i = 0; i < steps; ++i) {
        aphi.step(dt);
    }
    const auto aphiEnd = std::chrono::steady_clock::now();
    std::cout << "A-Phi pusher: " << std::chrono::duration_cast<std::chrono::microseconds>(aphiEnd - aphiBegin).count()*1e-6 <<"s \n";
#endif

#ifndef APHI_ONLY
    auto boris = BorisPusher(ef, mf);
    boris.setElectronInfo(offset+rLarmor*std::cos(startAngle), rLarmor*std::sin(startAngle),0, 0, u, 0);
    const auto borisBegin = std::chrono::steady_clock::now();
    for (int64_t i = 0; i < steps; ++i) {
        boris.step(dt);
    }
    const auto borisEnd = std::chrono::steady_clock::now();
    std::cout << "Boris pusher: " << std::chrono::duration_cast<std::chrono::microseconds>(borisEnd - borisBegin).count()*1e-6 <<"s \n";
#endif

    return 0;
}
