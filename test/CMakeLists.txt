
set(CMAKE_BUILD_TYPE Debug)

# The test files I care only (create a list)
set(testsources
	test3d_1.cpp
	# test3d_2.cpp
	# test3d_3.f90
	test3d_4.cpp
	test3d_5.cpp
	# test_MIC3D_8.cpp
	# test3d_6.c
	# test_get_elem_nodes.cpp
	test_ell_1.cpp
	test_ell_2.cpp
	# test_ell_mvp_openacc.cpp
	# test_cg.cpp
	# test_print_vtu_1.cpp
	# test_omp.cpp
	test_material.cpp
	test_damage.cpp
	test_util_1.cpp
	# test_A0.cpp
	# test_restart.cpp
	# benchmark-mic-1.cpp
	benchmark-mic-2.cpp
	# benchmark-mic-3.cpp
	# benchmark-mic-4.cpp
	benchmark-elastic.cpp
	benchmark-plastic.cpp
	benchmark-damage.cpp
	)

# Iterate over the list above
foreach (testfile ${testsources})

	# Delete File extensions (test_i.cpp -> test_i)
	get_filename_component(testname ${testfile} NAME_WE)

	# Create the executable test_i and link.
	add_executable(${testname} ${testfile})
	target_link_libraries(${testname} micropp)

	if(ENABLE_CUDA)
	  set_source_files_properties(${PROJECT_SOURCE_DIR}/test/${testfile} PROPERTIES LANGUAGE CUDA)
	endif()

endforeach ()

# As some tests will require commands and could be executed in loops we add the
# tests individually here.

add_test(NAME test3d_1 COMMAND test3d_1 5 0 10)
#add_test(NAME test3d_3 COMMAND test3d_3 5 5 5 10)
add_test(NAME test3d_4 COMMAND test3d_4 5 3 10)
add_test(NAME test3d_5 COMMAND test3d_5 5 5 5 2 10)
add_test(NAME test_ell_1 COMMAND test_ell_1)
add_test(NAME test_ell_2 COMMAND test_ell_2)
add_test(NAME test_util_1 COMMAND test_util_1)
add_test(NAME test_material COMMAND test_material 5)
add_test(NAME benchmark-elastic COMMAND benchmark-elastic)
add_test(NAME benchmark-plastic COMMAND benchmark-plastic)
add_test(NAME benchmark-damage COMMAND benchmark-damage)
add_test(NAME test_damage COMMAND test_damage 10)

#set_property(TARGET test3d_3 PROPERTY LINKER_LANGUAGE Fortran)
