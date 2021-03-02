cmake_minimum_required(VERSION 3.19)

# specify the compiler
set(CMAKE_C_COMPILER "gcc-10")
set(CMAKE_CXX_COMPILER "/usr/bin/g++-10")


# set the project name
project("TemplatizedRungekuttaCPP20")


#This is what brings rk45 into the project
include_directories(rk45)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++20")


# add the executable
add_executable(NewtonTester newton_test.cpp)

# test directory
add_subdirectory(tests)
