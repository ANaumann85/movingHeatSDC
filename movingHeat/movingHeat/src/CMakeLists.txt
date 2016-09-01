add_library(heat Heat.cc)

add_executable(heat_imexEuler heat_imexEuler.cc)
target_link_libraries(heat_imexEuler heat)
add_executable(heat_ros2 heat_ros2.cc)
target_link_libraries(heat_ros2 heat)
add_executable(printSrcPoint printSrcPoint.cc)
target_link_libraries(printSrcPoint heat)
add_executable(getEigenvalues getEigenvalues.cc)
target_link_libraries(getEigenvalues heat)

target_link_dune_default_libraries("heat_imexEuler")
target_link_dune_default_libraries("heat_ros2")
