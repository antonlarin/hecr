cmake_minimum_required(VERSION 3.1)
project(hecr)

macro(eval_bool_expression VAR)
    if(${ARGN})
        set(${VAR} 1)
    else()
        set(${VAR} 0)
    endif()
endmacro()

eval_bool_expression(USING_MSVC ${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
eval_bool_expression(USING_ICC ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
eval_bool_expression(USING_APPLECLANG
    ${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")

if(USING_MSVC)
    add_definitions(-D_USE_MATH_DEFINES)
endif()

if (USING_ICC)
    set(COMMON_COMPILE_OPTIONS "-std=c++11")
endif()

set(COMMON_DIR ${PROJECT_SOURCE_DIR}/common)
set(COMMON_SRC ${COMMON_DIR}/algorithms.cpp)

set(PARALLELISM_VARIANT seq CACHE STRING
    "select whether to build sequential, OpenMP, TBB or Cilk Plus version")

set(VALID_VARIANTS_REGEX "(seq|omp|tbb|cilk)")
if(NOT ${PARALLELISM_VARIANT} MATCHES ${VALID_VARIANTS_REGEX})
    message(FATAL_ERROR
        "\"${PARALLELISM_VARIANT}\" is not a valid variant value. \
        Parallelism variant can be either one of:
    seq (default)
    omp
    tbb
    cilk")
    return()
endif()

add_subdirectory(${PROJECT_SOURCE_DIR}/${PARALLELISM_VARIANT})
include_directories(${COMMON_DIR} ${VARIANT_INCLUDE_DIRS})

set(EXECUTABLE_NAME ${PROJECT_NAME}_${PARALLELISM_VARIANT})
add_executable(${EXECUTABLE_NAME} ${COMMON_DIR}/main.cpp ${COMMON_SRC}
    ${VARIANT_SRC})

if(NOT USING_ICC)
    target_compile_features(${EXECUTABLE_NAME} PRIVATE cxx_auto_type)
endif()
target_link_libraries(${EXECUTABLE_NAME} ${VARIANT_LIBS})
set_property(TARGET ${EXECUTABLE_NAME} PROPERTY COMPILE_OPTIONS
    ${COMMON_COMPILE_OPTIONS} ${VARIANT_COMPILE_OPTIONS})

# Testing support
set(GTEST_DIR ${PROJECT_SOURCE_DIR}/3rdparty/gtest)
add_subdirectory(${GTEST_DIR} ${CMAKE_BINARY_DIR}/gtest)
if(USING_MSVC)
    option(gtest_force_shared_crt ON)
elseif(USING_APPLECLANG)
    add_definitions(-DGTEST_USE_OWN_TR1_TUPLE=1)
endif()

set(TESTS_SRC ${COMMON_DIR}/tests.cpp)

include_directories(SYSTEM ${GTEST_DIR}/include ${GTEST_DIR})

set(TESTS_EXECUTABLE_NAME "${EXECUTABLE_NAME}_tests")
add_executable(${TESTS_EXECUTABLE_NAME} ${TESTS_SRC} ${COMMON_SRC}
    ${VARIANT_SRC})

if(NOT USING_ICC)
    target_compile_features(${EXECUTABLE_NAME} PRIVATE cxx_auto_type)
endif()
target_link_libraries(${TESTS_EXECUTABLE_NAME} gtest_main ${VARIANT_LIBS})
set_property(TARGET ${TESTS_EXECUTABLE_NAME} PROPERTY COMPILE_OPTIONS
    ${COMMON_COMPILE_OPTIONS} ${VARIANT_COMPILE_OPTIONS})

