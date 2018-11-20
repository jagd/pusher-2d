TEST_REGISTER(MagField
    test-magfield.cc
)
TEST_LINK_LIBRARY(MagField ${PROJECT_NAME})


TEST_REGISTER(Pushers3D
    test-pushers-3d.cc
)
TEST_LINK_LIBRARY(Pushers3D ${PROJECT_NAME})


TEST_REGISTER(APhi
    test-aphi.cc
)
TEST_LINK_LIBRARY(APhi ${PROJECT_NAME})
