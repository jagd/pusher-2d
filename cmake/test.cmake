enable_testing()
macro(TEST_REGISTER)
    set(_test_name ${ARGV0})
    string(TOLOWER "${_test_name}" _test_exec_tail)
    set(_test_exec "test_${_test_exec_tail}")
    list(APPEND unit_tests ${_test_exec})
    set(_list_var "${ARGN}")
    list(REMOVE_AT _list_var 0)
    add_test(NAME ${_test_name} COMMAND ${_test_exec} WORKING_DIRECTORY ${UNIT_TEST_BIN_OUTPUT_DIR})
    add_executable(${_test_exec} ${_list_var})
    target_compile_definitions(${_test_exec} PRIVATE
        BOOST_TEST_DYN_LINK
        BOOST_TEST_MODULE=${_test_name}
    )
    target_link_libraries(${_test_exec} ${GTEST_BOTH_LIBRARIES})
endmacro()

macro(TEST_LINK_LIBRARY)
    set(_test_name ${ARGV0})
    string(TOLOWER "${_test_name}" _test_exec_tail)
    set(_test_exec "test_${_test_exec_tail}")
    set(_list_var "${ARGN}")
    list(REMOVE_AT _list_var 0)
    target_link_libraries(${_test_exec} ${_list_var})
endmacro()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${unit_tests})
