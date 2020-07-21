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
message(STATUS "ROOT_INCLUDE_DIRS: " ${ROOT_INCLUDE_DIRS})

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include("${ROOT_USE_FILE}")

# Hack introduced to support build on Irinas laptop
if (NOT EXISTS ${ROOT_INCLUDE_DIRS})
  message(STATUS "ROOT_INCLUDE_DIR path does not exist - trying ROOTSYS")
  if (NOT EXISTS $ENV{ROOTSYS})
    message(FATAL_ERROR " ROOTSYS path doesn't exist. Unable to proceed.")
  endif()
  set(ROOT_INCLUDE_DIRS $ENV{ROOTSYS}/include)
  message(STATUS setting ROOT_INCLUDE_DIRS to ${ROOT_INCLUDE_DIRS})
endif()



#=============================================================================
# Include subdirectories
#=============================================================================
add_subdirectory(src)

#=============================================================================
# cmake summary
#=============================================================================
message(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})
message(STATUS "ROOT_INCLUDE_DIRS: " ${ROOT_INCLUDE_DIRS})
message(STATUS "ROOT_USE_FILE: " ${ROOT_USE_FILE})
message(STATUS "ROOT_SYS: " $ENV{ROOTSYS})
