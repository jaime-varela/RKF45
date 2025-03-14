cmake_minimum_required(VERSION 3.19)

# specify the compiler
find_program(GCC_PATH NAMES gcc-13 gcc-12 gcc-11 gcc-10 gcc)
if (GCC_PATH)
    set(CMAKE_C_COMPILER "${GCC_PATH}" CACHE STRING "" FORCE)
else()
    message(FATAL_ERROR "No suitable GCC compiler found (10+)")
endif()
# Specify the C++ compiler
find_program(GPP_PATH NAMES g++-13 g++-12 g++-11 g++-10 g++)
if (GPP_PATH)
    set(CMAKE_CXX_COMPILER "${GPP_PATH}" CACHE STRING "" FORCE)
else()
    message(FATAL_ERROR "No suitable G++ compiler found (10+)")
endif()


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
