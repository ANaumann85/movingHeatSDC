project(mrsdc)
cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 11)

add_executable(runsdc src/runsdc.cc)
add_executable(runsdc_inhom src/runsdc_inhom.cc)
add_executable(runsdc_inhom_mass src/runsdc_inhom_mass.cc)
add_executable(testCollocation src/testCollocation.cc)
add_executable(cmp_sdc_inhom src/cmp_sdc_inhom.cc)
add_executable(runsdc_model src/runsdc_model.cc)
add_executable(test_sdc_model src/test_sdc_model.cc)
add_executable(test_mrsdc_g0 src/test_mrsdc_g0.cc)

add_executable(runsdc_standard src/runsdc_standard.cc)
add_executable(stability_model_sdc src/stability_model_sdc.cc)
target_link_libraries(stability_model_sdc lapack)
add_executable(stability_model_sdc_standard src/stability_model_sdc_standard.cc)
target_link_libraries(stability_model_sdc_standard lapack)

enable_testing()
add_test(cmp_sdc_inhom cmp_sdc_inhom)
add_test(runsdc_inhom runsdc_inhom)
add_test(runsdc_inhom_mass runsdc_inhom_mass)
add_test(runsdc runsdc)
add_test(test_sdc_model test_sdc_model)
