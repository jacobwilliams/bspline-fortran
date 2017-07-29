#
# Basic CMake configuration file for bspline-fortran
#

cmake_minimum_required (VERSION 2.8)

set(PROJECT_DESCRIPTION "Multidimensional B-Spline Interpolation of Data on a Regular Grid")
set(PROJECT_URL "https://github.com/jacobwilliams/bspline-fortran")

project (bspline-fortran NONE)
enable_language (Fortran)

file(GLOB SOURCES src/*.f90)
add_library(${PROJECT_NAME} ${SOURCES})
