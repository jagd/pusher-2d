if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(CMAKE_CONFIGURATION_TYPES Debug Release RelWithDebInfo Coverage CACHE TYPE INTERNAL FORCE )
    set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
    find_program(GCOV gcov)
    find_program(LCOV lcov)
    find_program(GENHTML genhtml)
    add_custom_target(coverage
        COMMAND ${LCOV} -z -d .
        COMMAND ${CMAKE_CTEST_COMMAND}
        COMMAND ${LCOV} -d . -c -o coverage_all.txt
        COMMAND ${LCOV} -r coverage_all.txt 'test/*' '/usr/*' -o coverage_src.txt
        COMMAND ${GENHTML} --demangle-cpp -o coverage_report coverage_src.txt
        COMMAND ${CMAKE_COMMAND} -E remove coverage_all.txt coverage_src.txt
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        DEPENDS ${unit_tests}
    )
endif()
