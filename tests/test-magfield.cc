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
}

TEST(LinearBzField, Validate_k0) {
    const double b = 1.3;
    const auto m = LinearBzField(b, 0);
    const double r = 30;
    const double z = 100.0;
    ASSERT_DOUBLE_EQ(b, m.bz(z, r));
    ASSERT_DOUBLE_EQ(r*r*b/2, r*m.aTheta(z, r));
}

TEST(BiUniformMagField, SimpleCheck) {
    ASSERT_NO_THROW(BiUniformMagField(0, -1, 1));
    const double z0 = 0.1;
    const auto m = BiUniformMagField(z0, -2, 2);
    ASSERT_EQ(m.bz(1, 1), 2);
    ASSERT_EQ(m.bz(-1, 100), -2);
}
