cmake_minimum_required(VERSION 3.1)
project(hecr_seq)
add_executable(hecr_seq main.cpp algorithms.cpp)
target_compile_features(hecr_seq PRIVATE cxx_auto_type)

# Testing support
set(GTEST_DIR ../3rdparty/gtest)
add_subdirectory(${GTEST_DIR} ${CMAKE_BINARY_DIR}/gtest)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    option(gtest_force_shared_crt ON)
elseif(APPLE)
    add_definitions(-DGTEST_USE_OWN_TR1_TUPLE=1)
endif()

include_directories(SYSTEM ${GTEST_DIR}/include ${GTEST_DIR})
add_executable(hecr_seq_tests tests.cpp algorithms.cpp)
target_link_libraries(hecr_seq_tests gtest_main)

