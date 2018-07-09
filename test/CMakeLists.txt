
# The test files I care only (create a list)
set(testsources
  test2d_1.cpp
  test3d_1.cpp
  test3d_7.cpp
  test3d_8.cpp
  test3d_3.f90)

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

# Add tests with no arguments.
add_test(test3d_3 test3d_3)

# Add a test that requires arguments
add_test(NAME test2d_1 COMMAND test2d_1 5 5 10)
add_test(NAME test3d_1 COMMAND test3d_1 5 5 5 1)
add_test(NAME test3d_7 COMMAND test3d_7 5 5 5 10)
add_test(NAME test3d_8 COMMAND test3d_8 5 5 5 10)
