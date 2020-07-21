cmake_minimum_required(VERSION  3.9 FATAL_ERROR)
project(dream_endcap_analysis)


#=============================================================================
# Augment CMake with our custom scripts
#=============================================================================
set(EXTRA_MODULES_DIR ${CMAKE_CURRENT_LIST_DIR}/cmake)


#=============================================================================
# General configuration
#=============================================================================
include(${EXTRA_MODULES_DIR}/EnsureBuildType.cmake)


#=============================================================================
# ROOT configuration
#=============================================================================
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT

if (DEFINED ENV{ROOTSYS})
  message(STATUS "**** Build using ROOTSYS ****")
  set(CMAKE_CXX_STANDARD 11)
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
  set(ROOT_INCLUDE_DIRS $ENV{ROOTSYS}/include)
  include_directories(${ROOT_INCLUDE_DIRS})
  link_directories($ENV{ROOTSYS}/lib)
else()
  message(STATUS "**** Build using CMAKE ****")
  find_package(ROOT REQUIRED)
  include("${ROOT_USE_FILE}")
endif()

if (NOT EXISTS ${ROOT_INCLUDE_DIRS})
  message(FATAL_ERROR "**** ROOT_INCLUDE_DIR path does not exist. ****")
endif()


#=============================================================================
# Include subdirectories
#=============================================================================
add_subdirectory(src)


#=============================================================================
# cmake summary
#=============================================================================
message(STATUS "**************************************************************")
message(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})
message(STATUS "ROOT_INCLUDE_DIRS: " ${ROOT_INCLUDE_DIRS})
message(STATUS "ROOT_USE_FILE: " ${ROOT_USE_FILE})
message(STATUS "ROOT_SYS: " $ENV{ROOTSYS})
message(STATUS "**************************************************************")
