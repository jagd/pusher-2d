# Planar large orbit
add_executable(planar-largeorbit-boris-demo planar-largeorbit-boris.cc)
target_compile_definitions(planar-largeorbit-boris-demo PRIVATE DEMO)
target_link_libraries(planar-largeorbit-boris-demo aphi-pusher)

add_executable(planar-largeorbit-boris planar-largeorbit-boris.cc)
target_link_libraries(planar-largeorbit-boris aphi-pusher)

add_executable(planar-largeorbit-aphi-demo planar-largeorbit-aphi.cc)
target_compile_definitions(planar-largeorbit-aphi-demo PRIVATE DEMO)
target_link_libraries(planar-largeorbit-aphi-demo aphi-pusher)

add_executable(planar-largeorbit-aphi planar-largeorbit-aphi.cc)
target_link_libraries(planar-largeorbit-aphi aphi-pusher)


# Planar small orbit
add_executable(planar-smallorbit-aphi-demo planar-smallorbit-aphi.cc)
target_compile_definitions(planar-smallorbit-aphi-demo PRIVATE DEMO)
target_link_libraries(planar-smallorbit-aphi-demo aphi-pusher)

add_executable(planar-smallorbit-aphi planar-smallorbit-aphi.cc)
target_link_libraries(planar-smallorbit-aphi aphi-pusher)

add_executable(planar-smallorbit-boris-demo planar-smallorbit-boris.cc)
target_compile_definitions(planar-smallorbit-boris-demo PRIVATE DEMO)
target_link_libraries(planar-smallorbit-boris-demo aphi-pusher)

add_executable(planar-smallorbit-boris planar-smallorbit-boris.cc)
target_link_libraries(planar-smallorbit-boris aphi-pusher)


# Performance measurement based on small orbit benchmarks
add_executable(planar-perf planar-perf.cc)
target_link_libraries(planar-perf aphi-pusher)

add_executable(planar-perf-boris planar-perf.cc)
target_compile_definitions(planar-perf-boris PRIVATE BORIS_ONLY)
target_link_libraries(planar-perf-boris aphi-pusher)

add_executable(planar-perf-aphi planar-perf.cc)
target_compile_definitions(planar-perf-aphi PRIVATE APHI_ONLY)
target_link_libraries(planar-perf-aphi aphi-pusher)


# Constant Ez
add_executable(const-ez const-ez.cc)
target_link_libraries(const-ez aphi-pusher)

# Linear oscillator
add_executable(mirrored-ez-demo-boris mirrored-ez.cc)
target_compile_definitions(mirrored-ez-demo-boris PRIVATE DEMO_BORIS)
target_link_libraries(mirrored-ez-demo-boris aphi-pusher)

add_executable(mirrored-ez-demo-aphi mirrored-ez.cc)
target_compile_definitions(mirrored-ez-demo-aphi PRIVATE DEMO_APHI)
target_link_libraries(mirrored-ez-demo-aphi aphi-pusher)

add_executable(mirrored-ez mirrored-ez.cc)
target_link_libraries(mirrored-ez aphi-pusher)

add_executable(mirrored-ez-very-relativistic mirrored-ez.cc)
target_compile_definitions(mirrored-ez-very-relativistic PRIVATE VERY_RELATIVISTIC)
target_link_libraries(mirrored-ez-very-relativistic  aphi-pusher)

# ConstBzEr
add_executable(const-bzer-largeorbit-demo const-bzer-largeorbit.cc)
target_compile_definitions(const-bzer-largeorbit-demo PRIVATE DEMO)
target_link_libraries(const-bzer-largeorbit-demo aphi-pusher)
