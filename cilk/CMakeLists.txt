cmake_minimum_required(VERSION 3.1)
project(hecr_cilk)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(USING_MSVC TRUE)
else()
    set(USING_MSVC FALSE)
endif()

if (USING_MSVC)
    add_definitions(-D_USE_MATH_DEFINES)
endif()

add_executable(hecr_cilk main.cpp algorithms.cpp solver.cpp)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set_property(TARGET hecr_cilk APPEND
        PROPERTY COMPILE_OPTIONS -std=c++11)
else()
    set_property(TARGET hecr_cilk APPEND
        PROPERTY COMPILE_OPTIONS -fcilkplus)
    target_link_libraries(hecr_cilk cilkrts)
    target_compile_features(hecr_cilk PRIVATE cxx_auto_type)
endif()

# Testing support
set(GTEST_DIR ../3rdparty/gtest)
add_subdirectory(${GTEST_DIR} ${CMAKE_BINARY_DIR}/gtest)
if (USING_MSVC)
    option(gtest_force_shared_crt ON)
elseif(APPLE)
    add_definitions(-DGTEST_USE_OWN_TR1_TUPLE=1)
endif()

include_directories(SYSTEM ${GTEST_DIR}/include ${GTEST_DIR})
add_executable(hecr_cilk_tests tests.cpp algorithms.cpp solver.cpp)
set_property(TARGET hecr_cilk_tests APPEND
    PROPERTY COMPILE_OPTIONS -fcilkplus)
target_link_libraries(hecr_cilk_tests gtest_main cilkrts)

