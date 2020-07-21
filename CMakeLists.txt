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
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

# Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include("${ROOT_USE_FILE}")


#=============================================================================
# Include subdirectories
#=============================================================================
add_subdirectory(src)

#=============================================================================
# cmake summary
#=============================================================================
message(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})
message(STATUS "ROOT INCLUDE: " ${ROOT_INCLUDE_DIRS})
message(STATUS "ROOT_USE_FILE: " ${ROOT_USE_FILE})
message(STATUS "ROOT_SYS: " $ENV{ROOTSYS})
