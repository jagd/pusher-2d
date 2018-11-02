#include <pusher/pusher.h>
#include <iostream>
#include <chrono>

#if defined(APHI_ONLY) || defined(BORIS_ONLY) || defined(LEAPFROG_ONLY) || defined(RK4_ONLY) 
#define SINGLE_PUSHER
#endif

int main()
{
    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<ConstBzField>(0);

    const int64_t steps = 20000000;
    const double dt = 1e-15;

#if !defined(SINGLE_PUSHER) || defined(LEAPFROG_ONLY)
    auto lf = LeapFrogPusher(ef, mf);
    lf.setElectronInfo(1, 0, 0, 0, 0, 0);
    const auto lfBegin = std::chrono::steady_clock::now();
    for (int64_t i = 0; i < steps; ++i) {
        lf.step(dt);
    }
    const auto lfEnd = std::chrono::steady_clock::now();
    std::cout << "LeapFrog pusher: " << std::chrono::duration_cast<std::chrono::microseconds>(lfEnd - lfBegin).count()*1e-6 <<"s \n";
#endif

#if !defined(SINGLE_PUSHER) || defined(BORIS_ONLY)
    auto boris = BorisPusher(ef, mf);
    boris.setElectronInfo(1, 0, 0, 0, 0, 0);
    const auto borisBegin = std::chrono::steady_clock::now();
    for (int64_t i = 0; i < steps; ++i) {
        boris.step(dt);
    }
    const auto borisEnd = std::chrono::steady_clock::now();
    std::cout << "Boris pusher: " << std::chrono::duration_cast<std::chrono::microseconds>(borisEnd - borisBegin).count()*1e-6 <<"s \n";
#endif

#if !defined(SINGLE_PUSHER) || defined(APHI_ONLY)
    auto aphi = APhiPusher(ef, mf);
    aphi.setElectronInfo(0, 1, 0, 0, aphi.pTheta(0, 1, 0), 1);
    const auto aphiBegin = std::chrono::steady_clock::now();
    for (int64_t i = 0; i < steps; ++i) {
        aphi.step(dt);
    }
    const auto aphiEnd = std::chrono::steady_clock::now();
    std::cout << "A-Phi pusher: " << std::chrono::duration_cast<std::chrono::microseconds>(aphiEnd - aphiBegin).count()*1e-6 <<"s \n";
#endif

#if !defined(SINGLE_PUSHER) || defined(RK4_ONLY)
    auto rk = RK4Pusher(ef, mf);
    rk.setElectronInfo(1, 0, 0, 0, 0, 0);
    const auto rkBegin = std::chrono::steady_clock::now();
    for (int64_t i = 0; i < steps; ++i) {
        rk.step(dt);
    }
    const auto rkEnd = std::chrono::steady_clock::now();
    std::cout << "RK4 pusher: " << std::chrono::duration_cast<std::chrono::microseconds>(rkEnd - rkBegin).count()*1e-6 <<"s \n";
#endif
    return 0;
}
