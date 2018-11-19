#include <pusher/pusher.h>
#include <iostream>
#include <chrono>

class ZeroField : public IStaticMagField, public IStaticEField
{
    // Inherited via IStaticMagField
    virtual double aTheta(double z, double r) const override
    {
        return 0.0;
    }
    virtual double aTheta2z(double z, double r) const override
    {
        return 0.0;
    }
    virtual double aTheta2r(double z, double r) const override
    {
        return 0.0;
    }
    virtual double br(double z, double r) const override
    {
        return 0.0;
    }
    virtual double bz(double z, double r) const override
    {
        return 0.0;
    }

    // Inherited via IStaticEField
    virtual double ez(double z, double r) const override
    {
        return 0.0;
    }
    virtual double er(double z, double r) const override
    {
        return 0.0;
    }
};


#if defined(APHI_ONLY) || defined(BORIS_ONLY) || defined(LEAPFROG_ONLY) || defined(RK4_ONLY)
#define SINGLE_PUSHER
#endif

int main()
{
    const auto zf = std::make_shared<ZeroField>();
    const int64_t steps = 200000000;
    const double dt = 1e-15;

#if !defined(SINGLE_PUSHER) || defined(LEAPFROG_ONLY)
    auto lf = LeapFrogPusher(zf, zf);
    lf.setElectronInfo(1, 0, 0, 0, 0, 0);
    const auto lfBegin = std::chrono::steady_clock::now();
    for (int64_t i = 0; i < steps; ++i) {
        lf.step(dt);
    }
    const auto lfEnd = std::chrono::steady_clock::now();
    std::cout << "LeapFrog pusher: " << std::chrono::duration_cast<std::chrono::microseconds>(lfEnd - lfBegin).count()*1e-6 <<"s \n";
#endif

#if !defined(SINGLE_PUSHER) || defined(BORIS_ONLY)
    auto boris = BorisPusher(zf, zf);
    boris.setElectronInfo(1, 0, 0, 0, 0, 0);
    const auto borisBegin = std::chrono::steady_clock::now();
    for (int64_t i = 0; i < steps; ++i) {
        boris.step(dt);
    }
    const auto borisEnd = std::chrono::steady_clock::now();
    std::cout << "Boris pusher: " << std::chrono::duration_cast<std::chrono::microseconds>(borisEnd - borisBegin).count()*1e-6 <<"s \n";
#endif

#if !defined(SINGLE_PUSHER) || defined(APHI_ONLY)
    auto aphi = APhiPusher(zf, zf);
    aphi.setElectronInfo(0, 1, 0, 0, aphi.pTheta(0, 1, 0), 1);
    const auto aphiBegin = std::chrono::steady_clock::now();
    for (int64_t i = 0; i < steps; ++i) {
        aphi.step(dt);
    }
    const auto aphiEnd = std::chrono::steady_clock::now();
    std::cout << "A-Phi pusher: " << std::chrono::duration_cast<std::chrono::microseconds>(aphiEnd - aphiBegin).count()*1e-6 <<"s \n";
#endif

#if !defined(SINGLE_PUSHER) || defined(RK4_ONLY)
    auto rk = RK4Pusher(zf, zf);
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
