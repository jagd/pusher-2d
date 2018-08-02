# Global C++ Compiler Settings

################################# Enable C++11 #################################
if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU
OR ${CMAKE_CXX_COMPILER_ID} STREQUAL Clang
OR ${CMAKE_CXX_COMPILER_ID} STREQUAL Intel
)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11")
endif()

################################# Strict Warnings ##############################
if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU
OR ${CMAKE_CXX_COMPILER_ID} STREQUAL Clang
OR ${CMAKE_CXX_COMPILER_ID} STREQUAL Intel
)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fomit-frame-pointer")
    add_definitions(-DM_PI=3.14159265358979323846264338327950288)
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Zi /W3")
    set(CMAKE_CXX_FLAGS_RLEASE "${CMAKE_CXX_FLAGS_RELEASE} /W3")
    add_definitions(-D_USE_MATH_DEFINES)
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
endif()
