#include <gtest/gtest.h>
#include <pusher/magfield.h>
#include <pusher/efield.h>
#include <pusher/pusher.h>
#include <memory>

TEST(BorisPusher, Ctor) {
    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<ConstBzField>(0);
    BorisPusher(ef, mf);
}


TEST(BorisPusher, ZeroFieldsWithoutMotion) {
    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<ConstBzField>(0);
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


TEST(BorisPusher, ZeroFieldsWithMotion) {
    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<ConstBzField>(0);
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


TEST(BorisPusher, PlainLargeOrbit) {
    const double B = 1.0; // f ~ 28 GHz
    const double u = 1e8; // f ~ 26 GHz considering gamma
    const double r = -M0/Q0/B * u;
    const double gamma = u2gamma(u);
    const double v = u/gamma;
    const double omega = v/r;

    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<ConstBzField>(B);
    auto pusher = BorisPusher(ef, mf);

    const double dt = 1e-13;
    const double startAngle = dt * omega/2;
    pusher.setElectronInfo(r*std::cos(startAngle), r*std::sin(startAngle),0, 0, u, 0);
    for (int i = 0; i < 10; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        ASSERT_NEAR(r, std::sqrt(p.x*p.x+p.y*p.y), r*1e-6);
    }
}


TEST(BorisPusher, PlainSmallOrbit) {
    const double B = 1.0; // f ~ 28 GHz
    const double u = 1e8; // f ~ 26 GHz considering gamma
    const double rLarmor = -M0/Q0/B * u; // ~ 0.5 mm
    const double gamma = u2gamma(u);
    const double v = u/gamma;
    const double omega = v/rLarmor;

    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<ConstBzField>(B);
    auto pusher = BorisPusher(ef, mf);

    const double dt = 1e-13;
    const double startAngle = dt * omega/2;
    const double offset = 2*rLarmor;
    pusher.setElectronInfo(offset+rLarmor*std::cos(startAngle), rLarmor*std::sin(startAngle),0, 0, u, 0);
    for (int i = 0; i < 10; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        const double rSoll = std::sqrt(offset*offset +rLarmor*rLarmor + 2*offset*rLarmor*std::cos(startAngle+omega*dt*(i+1)));
        ASSERT_NEAR(rSoll, std::sqrt(p.x*p.x+p.y*p.y), rLarmor*1e-6);
    }
}


TEST(BorisPusher, MirroredEzField) {
    const auto ef = std::make_shared<MirroredEzField>(10e3/10e-2);
    const auto mf = std::make_shared<ConstBzField>(0);
    auto pusher = BorisPusher(ef, mf);

    const double zInit = -10e-2;
    const double totalEnergy = Q0*ef->pot(zInit, 0);
    const double refZ = -0.089871131212472837868; // for 10 kV / 10 cm
    const int pseudoOrder = 30;
    int64_t maxSteps = static_cast<int64_t>(1) << pseudoOrder;
    const double minDt = 1e-18;
    const double r = 1e-2;
    const int orderOffset = 15;
    int64_t steps = maxSteps >> orderOffset;
    double dt = minDt*(1 << orderOffset);
    pusher.setElectronInfo(r, 0, zInit, 0, 0, 0);
    for (int i = 0; i < steps; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        const double potEnergy = Q0*ef->pot(p.z, 0);
        const double kinEnergy = (pusher.gamma()-1.0)*(M0*C0*C0);
        ASSERT_NEAR(totalEnergy/Q0, (potEnergy+kinEnergy)/Q0, std::abs(totalEnergy/Q0)*1e-5);
    }
    ASSERT_NEAR(refZ, pusher.pos().z, std::abs(refZ*1e-5));
}


TEST(BorisPusher, ConstEzField) {
    const auto ef = std::make_shared<ConstEzField>(-10e3/10e-2);
    const auto mf = std::make_shared<ConstBzField>(0);
    auto pusher = BorisPusher(ef, mf);

    const double zInit = -10e-2;
    const double totalEnergy = Q0*ef->pot(zInit, 0);
    const double refZ = -0.089871131212472837868; // for 10 kV / 10 cm
    const int pseudoOrder = 30;
    int64_t maxSteps = static_cast<int64_t>(1) << pseudoOrder;
    const double minDt = 1e-18;
    const double r = 1e-2;
    const int orderOffset = 15;
    int64_t steps = maxSteps >> orderOffset;
    double dt = minDt*(1 << orderOffset);
    pusher.setElectronInfo(r, 0, zInit, 0, 0, 0);
    for (int i = 0; i < steps; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        const double potEnergy = Q0*ef->pot(p.z, 0);
        const double kinEnergy = (pusher.gamma()-1.0)*(M0*C0*C0);
        ASSERT_NEAR(totalEnergy/Q0, (potEnergy+kinEnergy)/Q0, std::abs(totalEnergy/Q0)*1e-5);
    }
    ASSERT_NEAR(refZ, pusher.pos().z, std::abs(refZ*1e-5));
}


TEST(BorisPusher, ConstErField) {
    const auto ef = std::make_shared<ConstErField>(-10e3/10e-2);
    const auto mf = std::make_shared<ConstBzField>(0);
    auto pusher = BorisPusher(ef, mf);

    const double rInit = 10e-2;
    const double totalEnergy = Q0*ef->pot(0, rInit);
    const double refR = rInit+0.010128868787527162132; // for 10 kV / 10 cm
    const int pseudoOrder = 30;
    int64_t maxSteps = static_cast<int64_t>(1) << pseudoOrder;
    const double minDt = 1e-18;
    const int orderOffset = 15;
    int64_t steps = maxSteps >> orderOffset;
    double dt = minDt*(1 << orderOffset);
    pusher.setElectronInfo(rInit, 0, 0, 0, 0, 0);
    for (int i = 0; i < steps; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        const double r = std::sqrt(p.x*p.x + p.y*p.y);
        const double potEnergy = Q0*ef->pot(0, r);
        const double kinEnergy = (pusher.gamma()-1.0)*(M0*C0*C0);
        ASSERT_NEAR(totalEnergy/Q0, (potEnergy+kinEnergy)/Q0, std::abs(totalEnergy/Q0)*1e-5);
    }
    const auto p = pusher.pos();
    ASSERT_NEAR(refR, std::sqrt(p.x*p.x + p.y*p.y), std::abs(refR*1e-5));
}
