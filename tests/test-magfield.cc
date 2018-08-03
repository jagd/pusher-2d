#include <gtest/gtest.h>
#include <pusher/magfield.h>

TEST(HomogeneousMagField, Ctor) {
    const auto m1 = HomogeneousMagField(0);
    const auto m2 = HomogeneousMagField(1.0);
    const auto m3 = HomogeneousMagField(-1.0);
}

