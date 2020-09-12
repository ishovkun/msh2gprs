# msh2gprs_root/tests/CMakeLists.txt

# shape function test polyhedral element direct with r0 subdivision
set(test_dir "${CMAKE_SOURCE_DIR}/examples/mech/pfem-direct-r0-hex")
add_test( pfem_direct_r0 "${PROJECT_BINARY_DIR}/msh2gprs" ${test_dir}/config.yaml)
add_test( pfem_direct_r0_compare ${CMAKE_COMMAND} -E compare_files
            ${test_dir}/output/gm_fem.txt ${test_dir}/gm_fem_test.txt)

# shape function test polyhedral element direct with gmsh subdivision
if(USE_GMSH)
  set(test_dir "${CMAKE_SOURCE_DIR}/examples/mech/pfem-direct-gmsh-hex")
  add_test( pfem_direct_gmsh "${PROJECT_BINARY_DIR}/msh2gprs" ${test_dir}/config.yaml)
  add_test( pfem_direct_gmsh_compare ${CMAKE_COMMAND} -E compare_files
            ${test_dir}/output/gm_fem.txt ${test_dir}/gm_fem_test.txt)
endif()

# beam bending : pebi grids + bc with location and value functions
set(test_dir "${CMAKE_SOURCE_DIR}/examples/mech/beam-shear-pebi")
add_test( beam_shear_pebi "${PROJECT_BINARY_DIR}/msh2gprs" ${test_dir}/config.yaml)
add_test( beam_shear_pebi_compare_bc ${CMAKE_COMMAND} -E compare_files
            ${test_dir}/output/gm_bcond.txt ${test_dir}/gm_bcond_test.txt)
add_test( beam_shear_pebi_compare_fem ${CMAKE_COMMAND} -E compare_files
            ${test_dir}/output/gm_fem.txt ${test_dir}/gm_fem_test.txt)

# pure mech, single dfm fracture + node BC
set(test_dir "${CMAKE_SOURCE_DIR}/examples/mech/dfm-shear")
add_test( dfm_shear_mech "${PROJECT_BINARY_DIR}/msh2gprs" ${test_dir}/config.yaml)
add_test( dfm_shear_mech_compare_bc ${CMAKE_COMMAND} -E compare_files
            ${test_dir}/output/gm_bcond.txt ${test_dir}/gm_bcond_test.txt)
add_test( dfm_shear_mech_compare_frac ${CMAKE_COMMAND} -E compare_files
            ${test_dir}/output/gm_DFM.txt ${test_dir}/gm_DFM_test.txt)

# pure flow, single edfm + simple well forced to connect
set(test_dir "${CMAKE_SOURCE_DIR}/examples/flow/simple-well-connect-edfm-2d")
add_test( well_edfm_force_connect "${PROJECT_BINARY_DIR}/msh2gprs" ${test_dir}/config.yaml)
add_test( well_edfm_force_connect_compare ${CMAKE_COMMAND} -E compare_files
            ${test_dir}/output/wells.txt ${test_dir}/wells_reference.txt)

# cartesian grid builder
set(test_dir "${CMAKE_SOURCE_DIR}/examples/flow/cartesian_grid_3d")
add_test( cartesian_grid_builder "${PROJECT_BINARY_DIR}/msh2gprs" ${test_dir}/config.yaml)
add_test( cartesian_grid_builder_compare ${CMAKE_COMMAND} -E compare_files
            ${test_dir}/output/flow_reservoir_mesh.vtk ${test_dir}/flow_reservoir_mesh_reference.vtk)
