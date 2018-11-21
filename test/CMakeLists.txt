
set(CMAKE_BUILD_TYPE Debug)

# The test files I care only (create a list)
set(testsources
	test2d_1.cpp
	test3d_1.cpp
	test3d_2.cpp
	test3d_3.f90
	test3d_4.cpp
	test3d_5.cpp
	test3d_6.c
	test3d_perf_1.cpp
	test_mictype_3.cpp
	test_MIC_QUAD_FIB_XZ.cpp
	test_get_elem_nodes.cpp
	test_ell_1.cpp
	test_ell_2.cpp
	test_ell_3.cpp
	conv3d_2lay.cpp
        test_print_vtu_1.cpp
	test_util_1.cpp)

# Iterate over the list above
foreach (testfile ${testsources})
	# Delete File extensions (test_i.cpp -> test_i)
	get_filename_component(testname ${testfile} NAME_WE)

	# Create the executable test_i and link.
	add_executable(${testname} ${testfile})
	target_link_libraries(${testname} micropp)
endforeach ()

# As some tests will require commands and could be executed in loops we add the
# tests individually here.

# Add a test that requires arguments
add_test(NAME test2d_1 COMMAND test2d_1 5 5 10)
add_test(NAME test3d_1 COMMAND test3d_1 5 5 5 1)
add_test(NAME test3d_2 COMMAND test3d_2 5 5 5 10)
add_test(NAME test3d_3 COMMAND test3d_3 5 5 5 10)
add_test(NAME test3d_4 COMMAND test3d_4 5 5 5 3 10)
add_test(NAME test3d_5 COMMAND test3d_5 5 5 5 2 10)
add_test(NAME test_ell_1 COMMAND test_ell_1)
add_test(NAME test_ell_2 COMMAND test_ell_2)
add_test(NAME test_ell_3 COMMAND test_ell_3 5 5)
add_test(NAME test_util_1 COMMAND test_util_1)
#add_test(NAME test_get_elem_nodes COMMAND test_get_elem_nodes)
