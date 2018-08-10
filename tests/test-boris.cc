#include <gtest/gtest.h>
#include <pusher/magfield.h>
#include <pusher/efield.h>
#include <pusher/pusher.h>
#include <memory>

TEST(Boris, Ctor) {
    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<HomogeneousMagField>(0);
    BorisPusher(ef, mf);
}


TEST(Boris, ZeroFieldsWithoutMotion) {
    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<HomogeneousMagField>(0);
    auto pusher = BorisPusher(ef, mf);
    pusher.setElectronInfo(0,0,0,0,0,0);
    for (int i = 0; i < 10; ++i) {
        pusher.step(1e-6);
        const auto p = pusher.pos();
        ASSERT_EQ(0, p.x);
        ASSERT_EQ(0, p.y);
        ASSERT_EQ(0, p.z);
    }
    pusher.setElectronInfo(1.0,0,0,0,0,0);
    for (int i = 0; i < 10; ++i) {
        pusher.step(1e-6);
        const auto p = pusher.pos();
        ASSERT_EQ(1.0, p.x);
        ASSERT_EQ(0, p.y);
        ASSERT_EQ(0, p.z);
    }
}


TEST(Boris, ZeroFieldsWithMotion) {
    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<HomogeneousMagField>(0);
    auto pusher = BorisPusher(ef, mf);
    const double dt = 1e-6;

    pusher.setElectronInfo(0,0,0,0,0,v2u(1.0));
    for (int i = 0; i < 10; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        ASSERT_EQ(0, p.x);
        ASSERT_EQ(0, p.y);
        ASSERT_DOUBLE_EQ(1.0*dt*(i+1), p.z);
    }

    pusher.setElectronInfo(0,1.0,0,0,0,v2u(1.0));
    for (int i = 0; i < 10; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        ASSERT_EQ(0, p.x);
        ASSERT_EQ(1.0, p.y);
        ASSERT_DOUBLE_EQ(1.0*dt*(i+1), p.z);
    }
}


TEST(Boris, PlainCyclotron) {
    const double B = 1.0; // f ~ 28 GHz
    const double u = 1e8; // f ~ 26 GHz considering gamma
    const double r = -M0/Q0/B * u;
    const double gamma = u2gamma(u);
    const double v = u/gamma;
    const double omega = v/r;

    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<HomogeneousMagField>(B);
    auto pusher = BorisPusher(ef, mf);

    double dt = 1e-13;
    double startAngle = dt * omega/2;
    pusher.setElectronInfo(r*std::cos(startAngle), r*std::sin(startAngle),0, 0, u, 0);
    for (int i = 0; i < 10; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        ASSERT_NEAR(r, std::sqrt(p.x*p.x+p.y*p.y), r*1e-6);
    }
}
