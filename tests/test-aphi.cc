#include <gtest/gtest.h>
#include <pusher/pusher.h>
#include <memory>

TEST(APhi, Ctor) {
    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<HomogeneousMagField>(0);
    APhiPusher(ef, mf);
}


TEST(APhi, ZeroFieldsWithoutMotion) {
    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<HomogeneousMagField>(0);
    auto pusher = APhiPusher(ef, mf);
    pusher.setElectronInfo(0,1,0,0,pusher.pTheta(0, 1, 0));
    for (int i = 0; i < 10; ++i) {
        pusher.step(1e-6);
        const auto p = pusher.pos();
        ASSERT_EQ(0, p.z);
        ASSERT_EQ(1, p.r);
    }
}


TEST(APhi, ZeroFieldsWithMotion) {
    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<HomogeneousMagField>(0);
    auto pusher = APhiPusher(ef, mf);
    const double dt = 1e-6;

    pusher.setElectronInfo(0,1,1.0,0,pusher.pTheta(0, 1, 0));
    for (int i = 0; i < 10; ++i) {
        pusher.step(1e-6);
        const auto p = pusher.pos();
        ASSERT_DOUBLE_EQ(1.0*dt*(i+1), p.z);
        ASSERT_EQ(1, p.r);
    }
}


TEST(APhi, PlainLargeOrbit) {
    const double B = 1.0; // f ~ 28 GHz
    const double u = 1e8; // f ~ 26 GHz considering gamma
    const double r = -M0/Q0/B * u;

    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<HomogeneousMagField>(B);
    auto pusher = APhiPusher(ef, mf);

    const double dt = 1e-13;
    pusher.setElectronInfo(0, r, 0, 0, pusher.pTheta(0, r, u));
    for (int i = 0; i < 10; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        ASSERT_EQ(0, p.z);
        ASSERT_NEAR(r, p.r, r*1e-12);
    }
}


TEST(APhi, PlainSmallOrbit) {
    const double B = 1.0; // f ~ 28 GHz
    const double u = 1e8; // f ~ 26 GHz considering gamma
    const double rLarmor = -M0/Q0/B * u;
    const double gamma = u2gamma(u);
    const double v = u/gamma;
    const double omega = v/rLarmor;

    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<HomogeneousMagField>(B);
    auto pusher = APhiPusher(ef, mf);

    const double offset = 2*rLarmor;
    const double dt = 1e-13;
    pusher.setElectronInfo(0, std::sqrt(offset*offset +rLarmor*rLarmor + 2*offset*rLarmor*std::cos(omega*dt*0.5)), 0, 0, pusher.pTheta(0, offset+rLarmor, u));
    for (int i = 0; i < 10; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        const double rSoll = std::sqrt(offset*offset +rLarmor*rLarmor + 2*offset*rLarmor*std::cos(omega*dt*(i+1.5)));
        ASSERT_NEAR(rSoll, p.r, rLarmor*1e-6);
    }
}
