# Math tests targets

include(FetchContent)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_Declare(
	googletest 
	GIT_REPOSITORY https://github.com/google/googletest.git 
	GIT_TAG main # Live at head
)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

get_property(existing_sources GLOBAL PROPERTY UNIT_TEST_SOURCES)

# Create a set for the Math unit tests targets
set(Math_UnitTest_TARGETS)
# Get a list of the .cpp files in the subdirectory for the unit tests
file(GLOB_RECURSE Math_UnitTest_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/unit/*.cpp")

list(APPEND existing_sources ${Math_UnitTest_Sources})

set_property(GLOBAL PROPERTY UNIT_TEST_SOURCES "${existing_sources}")

get_property(UNIT_TEST_TARGETS GLOBAL PROPERTY UNIT_TEST_TARGETS)

# Add each source file as a test target
foreach(TEST_SOURCE ${Math_UnitTest_SOURCES})
	get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
	add_executable(${TEST_NAME} EXCLUDE_FROM_ALL ${TEST_SOURCE})
	target_link_libraries(${TEST_NAME} PRIVATE GTest::gtest_main UtilitiesStatic)
	set_target_properties(${TEST_NAME} PROPERTIES INSTALLABLE OFF)
	list(APPEND Math_UnitTests_TARGETS ${TEST_NAME})
	list(APPEND UNIT_TEST_TARGETS ${TEST_NAME})
	gtest_discover_tests(${TEST_NAME} PROPERTIES)
endforeach()

add_custom_target(Math_UnitTests DEPENDS ${Math_UnitTests_TARGETS})
# Create an executable for the custom target, such that the IDEs can see it as a runnable target
add_executable(run_Math_UnitTests EXCLUDE_FROM_ALL ${Math_UnitTest_SOURCES})
# Link the executable with GTest and the TokenValueParser library
target_link_libraries(run_Math_UnitTests PRIVATE GTest::gtest_main UtilitiesStatic)
set_target_properties(run_Math_UnitTests PROPERTIES INSTALLABLE OFF)

set_property(GLOBAL PROPERTY UNIT_TEST_TARGETS "${UNIT_TEST_TARGETS}")