#include <gtest/gtest.h>
#include <pusher/pusher.h>
#include <memory>

TEST(APhiPusher, Ctor) {
    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<ConstBzField>(0);
    APhiPusher(ef, mf);
}


TEST(APhiPusher, ZeroFieldsWithoutMotion) {
    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<ConstBzField>(0);
    auto pusher = APhiPusher(ef, mf);
    pusher.setElectronInfo(0,1,0,0,pusher.pTheta(0, 1, 0), 1.0);
    for (int i = 0; i < 10; ++i) {
        pusher.step(1e-6);
        const auto p = pusher.pos();
        ASSERT_EQ(0, p.z);
        ASSERT_EQ(1, p.r);
    }
}


TEST(APhiPusher, ZeroFieldsWithMotion) {
    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<ConstBzField>(0);
    auto pusher = APhiPusher(ef, mf);
    const double dt = 1e-6;
    const double v = 1.0;
    const double gamma = 1 / std::sqrt(1 - v * v / C0 / C0);
    pusher.setElectronInfo(0, 1, v*gamma, 0, pusher.pTheta(0, 1, 0), gamma);
    for (int i = 0; i < 10; ++i) {
        pusher.step(1e-6);
        const auto p = pusher.pos();
        ASSERT_DOUBLE_EQ(1.0*dt*(i+1), p.z);
        ASSERT_EQ(1, p.r);
    }
}


TEST(APhiPusher, PlainLargeOrbit) {
    const double B = 1.0; // f ~ 28 GHz
    const double u = 1e8; // f ~ 26 GHz considering gamma
    const double r = -M0/Q0/B * u;

    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<ConstBzField>(B);
    auto pusher = APhiPusher(ef, mf);

    const double dt = 1e-13;
    pusher.setElectronInfo(0, r, 0, 0, pusher.pTheta(0, r, u), std::sqrt(1+u*u/C0/C0));
    for (int i = 0; i < 10; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        ASSERT_EQ(0, p.z);
        ASSERT_NEAR(r, p.r, r*1e-12);
    }
}


TEST(APhiPusher, PlainSmallOrbit) {
    const double B = 1.0; // f ~ 28 GHz
    const double u = 1e8; // f ~ 26 GHz considering gamma
    const double rLarmor = -M0/Q0/B * u;
    const double gamma = u2gamma(u);
    const double v = u/gamma;
    const double omega = v/rLarmor;

    const auto ef = std::make_shared<ZeroEField>();
    const auto mf = std::make_shared<ConstBzField>(B);
    auto pusher = APhiPusher(ef, mf);

    const double offset = 2*rLarmor;
    const double dt = 1e-13;
    pusher.setElectronInfo(
        0, std::sqrt(offset*offset +rLarmor*rLarmor + 2*offset*rLarmor*std::cos(omega*dt*0.5)),
        0, 0,
        pusher.pTheta(0, offset+rLarmor, u),
        std::sqrt(1+u*u/C0/C0)
    );
    for (int i = 0; i < 10; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        const double rSoll = std::sqrt(offset*offset +rLarmor*rLarmor + 2*offset*rLarmor*std::cos(omega*dt*(i+1.5)));
        ASSERT_NEAR(rSoll, p.r, rLarmor*1e-6);
    }
}


TEST(APhiPusher, MirroredEzField) {
    const auto ef = std::make_shared<MirroredEzField>(10e3/10e-2);
    const auto mf = std::make_shared<ConstBzField>(0);
    auto pusher = APhiPusher(ef, mf);

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
    pusher.setElectronInfo(zInit, r, 0, 0, pusher.pTheta(zInit, r, 0), 1.0);
    for (int i = 0; i < steps; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        const double potEnergy = Q0*ef->pot(p.z, 0);
        const double kinEnergy = (pusher.gammaCurrent()-1.0)*(M0*C0*C0);
        ASSERT_NEAR(
            totalEnergy/Q0,
            (potEnergy+kinEnergy)/Q0,
            std::abs(totalEnergy/Q0)*1e-12
        );
    }
    ASSERT_NEAR(refZ, pusher.pos().z, std::abs(refZ*1e-5));
}


TEST(APhiPusher, ConstEzField) {
    const auto ef = std::make_shared<ConstEzField>(-10e3/10e-2);
    const auto mf = std::make_shared<ConstBzField>(0);
    auto pusher = APhiPusher(ef, mf);

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
    pusher.setElectronInfo(zInit, r, 0, 0, pusher.pTheta(zInit, r, 0), 1.0);
    for (int i = 0; i < steps; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        const double potEnergy = Q0*ef->pot(p.z, 0);
        const double kinEnergy = (pusher.gammaCurrent()-1.0)*(M0*C0*C0);
        ASSERT_NEAR(
            totalEnergy/Q0,
            (potEnergy+kinEnergy)/Q0,
            std::abs(totalEnergy/Q0)*1e-12
        );
    }
    ASSERT_NEAR(refZ, pusher.pos().z, std::abs(refZ*1e-5));
}


TEST(APhiPusher, ConstErField) {
    const auto ef = std::make_shared<ConstErField>(-10e3/10e-2);
    const auto mf = std::make_shared<ConstBzField>(0);
    auto pusher = APhiPusher(ef, mf);

    const double rInit = 10e-2;
    const double totalEnergy = Q0*ef->pot(0, rInit);
    const double refR = rInit+0.010128868787527162132; // for 10 kV / 10 cm
    const int pseudoOrder = 30;
    int64_t maxSteps = static_cast<int64_t>(1) << pseudoOrder;
    const double minDt = 1e-18;
    const int orderOffset = 15;
    int64_t steps = maxSteps >> orderOffset;
    double dt = minDt*(1 << orderOffset);
    pusher.setElectronInfo(0, rInit, 0, 0, pusher.pTheta(0, rInit, 0), 1.0);
    for (int i = 0; i < steps; ++i) {
        pusher.step(dt);
        const auto p = pusher.pos();
        const double potEnergy = Q0*ef->pot(0, p.r);
        const double kinEnergy = (pusher.gammaCurrent()-1.0)*(M0*C0*C0);
        ASSERT_NEAR(
            totalEnergy/Q0,
            (potEnergy+kinEnergy)/Q0,
            std::abs(totalEnergy/Q0)*1e-12
        );
    }
    ASSERT_NEAR(refR, pusher.pos().r, std::abs(refR*1e-5));
}

static void auxConstBzEr(bool useDegradedLinearBzField = false)
{
	const double Bz = 1.0;
	const double Er = -10e3 / 10e-2;
	const double gamma = 1.2;
	const double r = 0.001131326056780128; // Mathematica
	const double v = gamma2v(gamma);
	const double omega = v / r;
	const double u = v * gamma;

	const auto ef = std::make_shared<ConstErField>(Er);
    const auto mf = useDegradedLinearBzField ?
		std::static_pointer_cast<IStaticMagField>(std::make_shared<LinearBzField>(Bz, 0)) :
		std::static_pointer_cast<IStaticMagField>(std::make_shared<ConstBzField>(Bz));
	auto aphi = APhiPusher(ef, mf);

	const double totalEnergy = Q0 * ef->pot(0, r) + (gamma - 1)*(M0*C0*C0);

	const double dt = 1e-15;
    const int64_t steps = static_cast<int>(M_PI/(omega*dt));
	aphi.setElectronInfo(0, r, 0, 0, aphi.pTheta(0, r, u), gamma);
	for (int i = 1; i <= steps; ++i) {
		aphi.step(dt);
	}
	const auto pAPhi = aphi.pos();
	const double kinEnergyAPhi = (aphi.gammaCurrent() - 1.0)*(M0*C0*C0);
	const double potEnergyAPhi = Q0 * ef->pot(pAPhi.z, pAPhi.r);
	ASSERT_NEAR(r, pAPhi.r, r*1e-8);
	ASSERT_NEAR(
		totalEnergy / Q0,
		(potEnergyAPhi + kinEnergyAPhi) / Q0,
		std::abs(totalEnergy / Q0 * 1e-8)
	);
}


TEST(APhiPusher, ConstBzErField)
{
	auxConstBzEr();
}


TEST(APhiPusher, LinearBzErField_k0)
{
	auxConstBzEr(true);
}