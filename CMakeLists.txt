cmake_minimum_required(VERSION 3.10)
project(LAB1)
add_library(LinLib source/linalg_realization.cpp)
target_include_directories(LinLib PUBLIC include)

