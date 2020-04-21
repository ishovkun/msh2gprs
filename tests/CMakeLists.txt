# msh2gprs_root/tests/CMakeLists.txt
# enable_testing()

# shape function test polyhedral element direct with r0 subdivision
set(test_dir "${CMAKE_SOURCE_DIR}/examples/mech/pfem-direct-r0-hex")
add_test( pfem_direct_r0 "${PROJECT_BINARY_DIR}/msh2gprs" ${test_dir}/config.yaml)
add_test( pfem_direct_r0_compare ${CMAKE_COMMAND} -E compare_files
            ${test_dir}/output/gm_fem.txt ${test_dir}/gm_fem_test.txt)

# shape function test polyhedral element direct with gmsh subdivision
set(test_dir "${CMAKE_SOURCE_DIR}/examples/mech/pfem-direct-gmsh-hex")
add_test( pfem_direct_gmsh "${PROJECT_BINARY_DIR}/msh2gprs" ${test_dir}/config.yaml)
add_test( pfem_direct_gmsh_compare ${CMAKE_COMMAND} -E compare_files
            ${test_dir}/output/gm_fem.txt ${test_dir}/gm_fem_test.txt)
