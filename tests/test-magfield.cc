#include <gtest/gtest.h>
#include <pusher/magfield.h>

TEST(ConstBzField, Ctor) {
    const auto m1 = ConstBzField(0);
    const auto m2 = ConstBzField(1.0);
    const auto m3 = ConstBzField(-1.0);
}

TEST(ConstBzField, Validate) {
    const double b = 1.3;
    const auto m = ConstBzField(b);
    const double r = 30;
    const double z = 100.0;
    ASSERT_DOUBLE_EQ(b, m.bz(z, r));
    ASSERT_DOUBLE_EQ(r*r*b/2, r*m.aTheta(z, r));
    double eps = 1e-6;
    ASSERT_NEAR(
        m.aTheta2r(z, r),
        (m.aTheta(z, r+eps/2) - m.aTheta(z, r-eps/2))/eps,
        m.aTheta2r(z, r)*1e-6
    );
}

TEST(LinearBzField, Validate_k0) {
    const double b = 1.3;
    const auto m = LinearBzField(b, 0);
    const double r = 30;
    const double z = 100.0;
    ASSERT_DOUBLE_EQ(b, m.bz(z, r));
    ASSERT_DOUBLE_EQ(r*r*b/2, r*m.aTheta(z, r));
    double eps = 1e-6;
    ASSERT_NEAR(
        m.aTheta2r(z, r),
        (m.aTheta(z, r+eps/2) - m.aTheta(z, r-eps/2))/eps,
        m.aTheta2r(z, r)*1e-6
    );
}

TEST(LinearBzField, Validate_alpha_Derivative) {
    const double b = 1.3;
    const auto m = LinearBzField(b, b);
    const double r = 30;
    const double z = 100.0;
    double eps = 1e-6;
    ASSERT_NEAR(
        m.aTheta2r(z, r),
        (m.aTheta(z, r+eps/2) - m.aTheta(z, r-eps/2))/eps,
        m.aTheta2r(z, r)*1e-6
    );
    ASSERT_NEAR(
        m.aTheta2z(z, r),
        (m.aTheta(z+eps/2, r) - m.aTheta(z-eps/2, r))/eps,
        m.aTheta2z(z, r)*1e-6
    );
}

TEST(BiUniformMagField, SimpleCheck) {
    ASSERT_NO_THROW(BiUniformMagField(0, -1, 1));
    const double z0 = 0.1;
    const auto m = BiUniformMagField(z0, -2, 2);
    ASSERT_EQ(m.bz(1, 1), 2);
    ASSERT_EQ(m.bz(-1, 100), -2);
}

TEST(BiUniformMagField, Validate_alpha_Derivative) {
    const auto m = BiUniformMagField(1.0, -1, 1);
    const double r = 30;
    const double z = 0.3;
    double eps = 1e-6;
    ASSERT_NEAR(
        m.aTheta2r(z, r),
        (m.aTheta(z, r+eps/2) - m.aTheta(z, r-eps/2))/eps,
        m.aTheta2r(z, r)*1e-6
    );
    ASSERT_NEAR(
        m.aTheta2z(z, r),
        (m.aTheta(z+eps/2, r) - m.aTheta(z-eps/2, r))/eps,
        m.aTheta2z(z, r)*1e-6
    );
}
