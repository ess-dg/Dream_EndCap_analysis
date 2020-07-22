include_directories(.)


set(read_bin_name "read_bin_irina")
set(analysis_name "analysis")
set(rootfile_name "CreateCDTRootFile")
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )

#
# Build read_bin
set(read_bin_src read_bin_irina.cpp)
add_executable( ${read_bin_name} ${read_bin_src} )


#
# analysis
set(analysis_src analysis_cdt.cpp)
add_executable(${analysis_name} ${analysis_src})
target_link_libraries(${analysis_name} RIO Core Tree Gpad)


#
# create tree
set(rootfile_src CreateCDTTree.cpp)
add_executable(${rootfile_name} ${rootfile_src})
target_link_libraries(${rootfile_name} RIO Core Tree)

install(TARGETS ${read_bin_name} ${analysis_name} ${rootfile_name} DESTINATION bin)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/endcap_lookup.root
        DESTINATION ${PROJECT_BINARY_DIR}/endcap_lookup.root)
