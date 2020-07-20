include_directories(.)

set(cdt_library)
add_library(${cdt_library} SHARED analysis_cdt.C CreateCDTTree.C)
target_link_libraries(${cdt_library} ${ROOT_LIBRARIES})

set(dea_src
   read_bin_irina.cpp
)

set(exec_name "read_bin_irina")
add_executable(
  ${exec_name} ${dea_src}
)

set_target_properties(${exec_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

install(TARGETS ${exec_name} DESTINATION bin)
