
# The test files I care only
set(testsources
  test_7.cpp
  test_8.cpp
  test_3.f90)

foreach (testfile ${testsources})
  # Parsing path and filename
  get_filename_component(testname ${testfile} NAME_WE)  # delete extension

  # Create the executable and link
  add_executable(${testname} ${testfile})
  target_link_libraries(${testname} microcpp)

  # Add a test with no arguments
  add_test(${testname} ${testname})
endforeach ()
