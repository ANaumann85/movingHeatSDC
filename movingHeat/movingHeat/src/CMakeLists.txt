add_library(heat Heat.cc)

add_executable(heat_imexEuler heat_imexEuler.cc)
target_link_libraries(heat_imexEuler heat)

target_link_dune_default_libraries("heat_imexEuler")
