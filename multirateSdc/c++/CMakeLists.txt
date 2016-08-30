project(mrsdc)
cmake_minimum_required(VERSION 3.5)

find_package(MTL REQUIRED)
include_directories(${MTL_INCLUDE_DIRS})
add_definitions(${MTL_CXX_DEFINITIONS})
add_executable(runsdc src/runsdc.cc)
