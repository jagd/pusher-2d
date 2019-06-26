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
add_executable(perf-zero perf-zero.cc)
target_link_libraries(perf-zero aphi-pusher)

add_executable(perf-zero-boris perf-zero.cc)
target_compile_definitions(perf-zero-boris PRIVATE BORIS_ONLY)
target_link_libraries(perf-zero-boris aphi-pusher)

add_executable(perf-zero-aphi perf-zero.cc)
target_compile_definitions(perf-zero-aphi PRIVATE APHI_ONLY)
target_link_libraries(perf-zero-aphi aphi-pusher)

add_executable(perf-zero-leapfrog perf-zero.cc)
target_compile_definitions(perf-zero-leapfrog PRIVATE LEAPFROG_ONLY)
target_link_libraries(perf-zero-leapfrog aphi-pusher)

add_executable(perf-zero-rk4 perf-zero.cc)
target_compile_definitions(perf-zero-rk4 PRIVATE RK4_ONLY)
target_link_libraries(perf-zero-rk4 aphi-pusher)

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
add_executable(const-bzer-largeorbit const-bzer-largeorbit.cc)
target_link_libraries(const-bzer-largeorbit aphi-pusher)

# ConstBzEz
add_executable(const-bzez-largeorbit-demo const-bzez-largeorbit.cc)
target_compile_definitions(const-bzez-largeorbit-demo PRIVATE DEMO)
target_link_libraries(const-bzez-largeorbit-demo aphi-pusher)

# Mixed field reflection demo
add_executable(mixed-reflection-demo mixed-reflection.cc)
target_compile_definitions(mixed-reflection-demo PRIVATE DEMO)
target_link_libraries(mixed-reflection-demo aphi-pusher)

# Constant Ez linear motion (Relativistic)
add_executable(ezbr-demo ezbr.cc)
target_compile_definitions(ezbr-demo PRIVATE DEMO)
target_link_libraries(ezbr-demo aphi-pusher)

add_executable(ezbr ezbr.cc)
target_link_libraries(ezbr aphi-pusher)

# cusp with pure magnetic field
add_executable(cusp-demo cusp.cc)
target_compile_definitions(cusp-demo PRIVATE DEMO)
target_link_libraries(cusp-demo aphi-pusher)

add_executable(cusp cusp.cc)
target_link_libraries(cusp aphi-pusher)
