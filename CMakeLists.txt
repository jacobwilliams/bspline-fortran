#
# Basic CMake configuration file for bspline-fortran
#

cmake_minimum_required(VERSION 3.11)

set( CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/install"
  CACHE STRING "Select where to install the library." )

project(bspline-fortran NONE)

# Define build options
option(BSPLINE_BUILD_TESTING "Build the test programs" ON)
option(BSPLINE_BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BSPLINE_BUILD_STATIC_LIBS "Build static libraries" ON)

# Real and Integer kinds
SET(BSPLINE_REAL_KIND "REAL64" CACHE STRING "Real kind parameter")
SET_PROPERTY(CACHE BSPLINE_REAL_KIND PROPERTY STRINGS REAL32 REAL64 REAL128)
if(${BSPLINE_REAL_KIND} MATCHES "REAL32")
  add_definitions(-DREAL32)
elseif(${BSPLINE_REAL_KIND} MATCHES "REAL64")
  add_definitions(-DREAL64)
elseif(${BSPLINE_REAL_KIND} MATCHES "REAL128")
  add_definitions(-DREAL128)
endif()
SET(BSPLINE_INT_KIND "INT32" CACHE STRING "Integer kind parameter")
SET_PROPERTY(CACHE BSPLINE_INT_KIND PROPERTY STRINGS INT8 INT16 INT32 INT64)
if(${BSPLINE_INT_KIND} MATCHES "INT8")
  add_definitions(-DINT8)
elseif(${BSPLINE_INT_KIND} MATCHES "INT16")
  add_definitions(-DINT16)
elseif(${BSPLINE_INT_KIND} MATCHES "INT32")
  add_definitions(-DINT32)
elseif(${BSPLINE_INT_KIND} MATCHES "INT64")
  add_definitions(-DINT64)
endif()

# Use MSVS folders to organize projects on windows
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set( LIB_NAME ${PROJECT_NAME} )
set( PROJECT_DESCRIPTION
  "Multidimensional B-Spline Interpolation of Data on a Regular Grid" )
set( PROJECT_URL "https://github.com/jacobwilliams/bspline-fortran" )
set( CMAKE_CONFIGURATION_TYPES "Debug" "Release" "MinSizeRel" "RelWithDebInfo" )
set( CMAKE_BUILD_TYPE "Release"
  CACHE STRING "Select which configuration to build." )

enable_language(Fortran)

set( SOURCES
  src/bspline_blas_module.F90
  src/bspline_defc_module.F90
  src/bspline_kinds_module.F90
  src/bspline_module.f90
  src/bspline_oo_module.f90
  src/bspline_sub_module.f90 )

set(INSTALL_MOD_DIR "lib")
set(INSTALL_LIB_DIR "${INSTALL_MOD_DIR}")
set(ABS_LIB_INSTALL_DIR "\${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}")
set(MODULE_DIR "${CMAKE_BINARY_DIR}/lib")

if( "${CMAKE_SYSTEM_NAME}" MATCHES "Darwin" )
  set( ENABLE_DYLIBS_USE_RPATH TRUE CACHE BOOL
    "Enable @rpath install name for dylibs" )
  mark_as_advanced( ENABLE_DYLIBS_USE_RPATH )
endif()

if( ENABLE_DYLIBS_USE_RPATH )
  set( CMAKE_MACOSX_RPATH TRUE )
else()
  set( CMAKE_INSTALL_NAME_DIR
    "${ABS_LIB_INSTALL_DIR}" )
endif()

if(BSPLINE_BUILD_SHARED_LIBS)
  add_library( ${LIB_NAME}        SHARED ${SOURCES} )

  set_target_properties(${LIB_NAME}
    PROPERTIES
    OUTPUT_NAME "${LIB_NAME}"
    if(NOT MSVC_IDE)
      PREFIX lib
    endif()
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

  # install library:
  install( TARGETS ${LIB_NAME}
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" )
  install( DIRECTORY "${MODULE_DIR}/" DESTINATION "${INSTALL_MOD_DIR}" )
endif()

if(BSPLINE_BUILD_STATIC_LIBS)
  add_library( ${LIB_NAME}-static STATIC ${SOURCES} )

  set_target_properties(${LIB_NAME}-static
    PROPERTIES
    OUTPUT_NAME "${LIB_NAME}"
    if(NOT MSVC_IDE)
      PREFIX lib
    endif()
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    Fortran_MODULE_DIRECTORY ${MODULE_DIR})

  # install library:
  install( TARGETS ${LIB_NAME}-static
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" )
  install( DIRECTORY "${MODULE_DIR}/" DESTINATION "${INSTALL_MOD_DIR}" )
endif()

# Windows settings:
if(MSVC_IDE)
   INCLUDE_DIRECTORIES("src")
   SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp")
endif()

if(BSPLINE_BUILD_TESTING)
  # Add ctest commands, to build all test programs (that do not depend on pyplot_module)
  # Thus, tests can all be run with `make test`
  enable_testing()
  foreach(prog_name
          # bspline_extrap_test
          bspline_test_2
          bspline_test
          # dbint4_test
          # knot_tests
          # speed_test
          # speed_test_oo
          stack_size_test_oo
          test_integrate
          test_oo
          test_regrid)
      add_executable(${prog_name} test/${prog_name}.f90)
      if(BSPLINE_BUILD_STATIC_LIBS)
        target_link_libraries(${prog_name} ${LIB_NAME}-static)
        target_include_directories(${prog_name} PRIVATE ${MODULE_DIR})
      elseif(BSPLINE_BUILD_SHARED_LIBS)
        target_link_libraries(${prog_name} ${LIB_NAME})
      else()
        message(FATAL_ERROR "No libraries to link the test programs to. Either shared or static libraries must be built.")
      endif()
      add_test(NAME ${prog_name}
              COMMAND ${prog_name})
  endforeach()
endif()
