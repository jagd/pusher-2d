#include <gtest/gtest.h>
#include <pusher/magfield.h>

TEST(HomogeneousMagField, Ctor) {
    const auto m1 = HomogeneousMagField(0);
    const auto m2 = HomogeneousMagField(1.0);
    const auto m3 = HomogeneousMagField(-1.0);
}

TEST(HomogeneousMagField, Validate) {
    const double b = 1.3;
    const auto m = HomogeneousMagField(b);
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
