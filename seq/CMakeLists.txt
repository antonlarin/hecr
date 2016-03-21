cmake_minimum_required(VERSION 3.1)
project(hecr_seq)
add_executable(hecr_seq main.cpp)
target_compile_features(hecr_seq PRIVATE cxx_auto_type)
