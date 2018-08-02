cmake_minimum_required(VERSION 2.8)

project(aphi-pusher)
add_executable(
    ${PROJECT_NAME}
    "src/main.cc"
    "src/magfield.h"
    "src/magfield.cc"
    )
