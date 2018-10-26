# Planar large orbit
add_executable(planar-largeorbit-demo planar-largeorbit.cc)
target_compile_definitions(planar-largeorbit-demo PRIVATE DEMO)
target_link_libraries(planar-largeorbit-demo aphi-pusher)

add_executable(planar-largeorbit planar-largeorbit.cc)
target_link_libraries(planar-largeorbit aphi-pusher)

# Planar small orbit
add_executable(planar-smallorbit-demo planar-smallorbit.cc)
target_compile_definitions(planar-smallorbit-demo PRIVATE DEMO)
target_link_libraries(planar-smallorbit-demo aphi-pusher)

add_executable(planar-smallorbit planar-smallorbit.cc)
target_link_libraries(planar-smallorbit aphi-pusher)

add_executable(planar-smallorbit-perioderror-demo planar-smallorbit-perioderror.cc)
target_link_libraries(planar-smallorbit-perioderror-demo aphi-pusher)

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

# ConstBzEz
add_executable(const-bzez-largeorbit-demo const-bzez-largeorbit.cc)
target_compile_definitions(const-bzez-largeorbit-demo PRIVATE DEMO)
target_link_libraries(const-bzez-largeorbit-demo aphi-pusher)

# Non-adiabatic: zero mag field
add_executable(nonadiabatic-cusp-demo nonadiabatic-cusp.cc)
target_compile_definitions(nonadiabatic-cusp-demo PRIVATE DEMO)
target_link_libraries(nonadiabatic-cusp-demo aphi-pusher)

# Linear moving (Relativistic)
add_executable(linear-relativistic linear-relativistic.cc)
target_compile_definitions(linear-relativistic PRIVATE DEMO)
target_link_libraries(linear-relativistic aphi-pusher)
