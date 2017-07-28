cmake_minimum_required (VERSION 2.8)
project (bspline-fortran NONE)
enable_language (Fortran)

file(GLOB SOURCES src/*.f90)
add_library(${PROJECT_NAME} ${SOURCES})
