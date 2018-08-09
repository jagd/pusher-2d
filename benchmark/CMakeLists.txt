add_executable(planar-largeorbit-boris-demo planar-largeorbit-boris.cc)
target_compile_definitions(planar-largeorbit-boris-demo PRIVATE DEMO)
target_link_libraries(planar-largeorbit-boris-demo aphi-pusher)

add_executable(planar-largeorbit-boris planar-largeorbit-boris.cc)
target_link_libraries(planar-largeorbit-boris aphi-pusher)
