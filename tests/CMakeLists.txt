# msh2gprs_root/tests/CMakeLists.txt
# enable_testing()

# dfm vertex splitting test
set(test_dir "${CMAKE_PROJECT_DIR}/examples/mech/dfm-shear")
add_test( dfm_shear_mech "${PROJECT_BINARY_DIR}/msh2gprs" ${test_dir}/conf.yaml)
add_test( dfm_shear_mech_compare ${CMAKE_COMMAND} -E compare_files
            ${test_dir}/dfm.vtk ${test_dir}/dfm1.vtk)
