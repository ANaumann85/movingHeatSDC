project(mrsdc)
cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 11)

add_executable(runsdc src/runsdc.cc)
add_executable(runsdc_inhom src/runsdc_inhom.cc)
add_executable(testCollocation src/testCollocation.cc)
