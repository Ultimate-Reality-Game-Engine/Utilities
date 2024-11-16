# CMakeList.txt : UltReality::Utilities

# CMake configuration and settings ****************************************************************
#**************************************************************************************************

cmake_minimum_required (VERSION 3.25)

project(UTILITIES VERSION 1.0.0)

# Google tests require mimimum c++ version 14
# setting minimum c++ version to 20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

# Define custom configuration types for Debug, Release, and Release-Dev
set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo" CACHE STRING "List of supported build modes." FORCE)

# Define compiler definitions to create OS specific code segments
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	add_compile_definitions(TARGET_OS=Windows _WIN_TARGET)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
	add_compile_definitions(TARGET_OS=Linux _LINUX_TARGET)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin") #MacOS
	add_compile_definitions(TARGET_OS=MacOS _MACOS_TARGET)
endif()

# Make sure CMake exports its compile commands so things like IntelliSense can detect them in editor
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Create options that are dependent onthis project being top level
option(Utilities_INSTALL "Enable project install" ${PROJECT_IS_TOP_LEVEL})
option(Utilities_BUILD_TESTS "Build the test suit" ${PROJECT_IS_TOP_LEVEL})

include(GNUInstallDirs)

# End CMake configuration and settings ************************************************************
#**************************************************************************************************

# Generated source ********************************************************************************
#**************************************************************************************************

# Create CMake referecne to a directory for configuration templates
set(CONFIG_TEMPLATES_DIR ${CMAKE_SOURCE_DIR}/config_templates)

# Add preprocessor definitions for conditional compilation based on configuration
#add_compile_definitions($<$<CONFIG:Debug>:_DEBUG> $<$<CONFIG:RelWithDebInfo>:_RELEASE_DEV> $<$<CONFIG:Release>:_RELEASE>)

# Create CMake reference to a directory for CMake generated files
set(GENERATED_DIR ${CMAKE_SOURCE_DIR}/generated)
# Create CMake reference to generated directory sub-directories for headers and source files
set(GENERATED_INCLUDE_DIR ${GENERATED_DIR}/include)
set(GENERATED_SRC_DIR ${GENERATED_DIR}/src)

# Create the generated directory and sub-directory structure if it doesn't exist
file(MAKE_DIRECTORY ${GENERATED_INCLUDE_DIR})
file(MAKE_DIRECTORY ${GENERATED_SRC_DIR})

# Function to process the configuration templates and generate source
function(generate_from_template TEMPLATE_DIR OUTPUT_DIR FILE_TYPE)
	file(GLOB TEMPLATES "${TEMPLATE_DIR}/*.in")
	foreach(TEMPLATE ${TEMPLATES})
		get_filename_component(FILENAME ${TEMPLATE} NAME_WE)
		set(OUTPUT_FILE ${OUTPUT_DIR}/${FILENAME}.${FILE_TYPE})
		configure_file(${TEMPLATE} ${OUTPUT_FILE} @ONLY)
	endforeach()
endfunction()

# Generate headers and source files from templates
generate_from_template(${CONFIG_TEMPLATES_DIR}/include ${GENERATED_INCLUDE_DIR} h)
generate_from_template(${CONFIG_TEMPLATES_DIR}/src ${GENERATED_SRC_DIR} cpp)

# End Generated Source ****************************************************************************
#**************************************************************************************************

# Target Creation *********************************************************************************
#**************************************************************************************************

# Create list of all source files to be used in this target
file(GLOB_RECURSE UtilitiesShared_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/Shared/*/src/*.cpp")
file(GLOB_RECURSE Shared_DIRS LIST_DIRECTORIES true ${CMAKE_CURRENT_SOURCE_DIR}/Shared/*)

set(UtilitiesShared_HEADERS "")
foreach(dir ${Shared_DIRS})
	if(IS_DIRECTORY "${dir}/include")
		list(APPEND UtilitiesShared_HEADERS "${dir}/include")
	endif()
endforeach()
list(APPEND UtilitiesShared_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/Shared/include")

#message(STATUS ${UtilitiesShared_HEADERS})

if(NOT DEFINED UTILITIES_SHARED_BUILD)
	set(UTILITIES_SHARED_BUILD ON)
endif()

if(UTILITIES_SHARED_BUILD)
	# Create the UtilitiesShared target
	if(UtilitiesShared_SOURCE)
		add_library(UtilitiesShared SHARED ${UtilitiesShared_SOURCE})

		target_compile_definitions(UtilitiesShared PRIVATE LIBRARY_EXPORTS)

		# Set up include directories for the shared library target
		foreach(header_dir ${UtilitiesShared_HEADERS})
			target_include_directories(UtilitiesShared PUBLIC 
				$<BUILD_INTERFACE:${header_dir}> 
 				$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> # This is used when the library is installed
			)
		endforeach()

		# Link libraries to the target
		#target_link_libraries()

		set_target_properties(UtilitiesShared PROPERTIES 
			VERSION ${PROJECT_VERSION} 
			SOVERSION ${PROJECT_VERSION_MAJOR}
		)
	endif()
endif()

file(GLOB_RECURSE UtilitiesStatic_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/Static/*/src/*.cpp")
file(GLOB_RECURSE Static_DIRS LIST_DIRECTORIES true ${CMAKE_CURRENT_SOURCE_DIR}/Static/*)

set(UtilitiesStatic_HEADERS "")
foreach(dir ${Static_DIRS})
	if(IS_DIRECTORY "${dir}/include")
		list(APPEND UtilitiesStatic_HEADERS "${dir}/include")
	endif()
endforeach()

if(NOT DEFINED UTILITIES_STATIC_BUILD)
	set(UTILITIES_STATIC_BUILD ON)
endif()

if(UTILITIES_STATIC_BUILD)
	# Create the UtilitiesStatic target
	if(UtilitiesStatic_SOURCE)
		add_library(UtilitiesStatic STATIC ${UtilitiesStatic_SOURCE})

		# Set up include directories for the shared library target
		foreach(header_dir ${UtilitiesStatic_HEADERS})
			target_include_directories(UtilitiesStatic PUBLIC 
				$<BUILD_INTERFACE:${header_dir}> 
 				$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> # This is used when the library is installed
			)
		endforeach()

		# Link libraries to the target
		#target_link_libraries()

		set_target_properties(UtilitiesStatic PROPERTIES 
			VERSION ${PROJECT_VERSION} 
			SOVERSION ${PROJECT_VERSION_MAJOR}
		)
	endif()
endif()

# End Target Creation *****************************************************************************
#**************************************************************************************************

# Installation and Packing Configuration **********************************************************
#**************************************************************************************************
if(Utilities_INSTALL)
	#if(UTILITIES_STATIC_BUILD) 
	#	set(Utilities_TARGETS UtilitiesShared UtilitiesStatic)
	#else()
	#	set(Utilities_TARGETS UtilitiesShared)
	#endif()

	# Install the targets
	install(
		TARGETS UtilitiesShared UtilitiesStatic
		EXPORT Utilities_Targets 
		ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # Static libraries/import libraries (.lib files for .dll linking) 
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # Shared libraries (.so) 
		RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # .exe or .dll 
		PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # Headers/include directories marked as PUBLIC 
		PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # Headers/include directories marked as PRIVATE
	)

	# Create the targets CMake file which contains the above definitions
	install(
		EXPORT Utilities_Targets 
		FILE Utilities_Targets.cmake 
		NAMESPACE UltReality:: 
		DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Utilities/cmake
	)

	# Install the actual includes
	install(
		DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ 
		DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	)

	# Generate and install the package version config files
	include(CMakePackageConfigHelpers)
	write_basic_package_version_file(
		"Utilities_ConfigVersion.cmake" 
		VERSION ${PROJECT_VERSION} 
		COMPATIBILITY SameMajorVersion
	)
	configure_package_config_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/cmake_config/Utilities_Config.cmake.in" 
		"${CMAKE_CURRENT_BINARY_DIR}/Utilities_Config.cmake" 
		INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Utilities/cmake
	)

	# Install the CMake config files
	install(
		FILES "${CMAKE_CURRENT_BINARY_DIR}/Utilities_ConfigVersion.cmake" 
		"${CMAKE_CURRENT_BINARY_DIR}/Utilities_Config.cmake" 
		DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Utilities/cmake
	)

	# Define Package install paths
	set(INCLUDEDIR_FOR_PKG_CONFIG "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
	set(LIBDIR_PKG_CONFIG "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

	# Create and install the package config file
	configure_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/cmake_config/Utilities.pc.in" 
		"${CMAKE_CURRENT_BINARY_DIR}/Utilities.pc" @ONLY
	)

	# Install the package config file
	install(
		FILES "${CMAKE_CURRENT_BINARY_DIR}/Utilities.pc" 
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
	)

	# A version that is often used to denote a specific build of the software, including revisions, builds, or other metadata
	set(PACKAGE_VERSION_BUILD "${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}")

	set(PACKAGE_VERSION "${PROJECT_VERSION}-${PACKAGE_VERSION_BUILD}")

	set(CPACK_PACKAGE_DIRECTORY "${CMAKE_SOURCE_DIR}/out/package")

	set(CPACK_PACKAGE_NAME "Utilities")
	set(CPACK_PACKAGE_VERSION "${PACKAGE_VERSION}")

	set(CPACK_PACKAGE_VENDOR "Andrew Todd")
	set(CPACK_PACKAGE_CONTACT "andrewdanieltodd@gmail.com")
	include(CPack)
	message(STATUS "PACKAGE_VERSION is: ${PACKAGE_VERSION}")
	message(STATUS "PACKAGE_FILE_NAME is: ${CPACK_PACKAGE_FILE_NAME}")
endif()

# End Installation and Packing Configuration ******************************************************
#**************************************************************************************************

# Create Unit Test Groups *************************************************************************
#**************************************************************************************************
if (Utilities_BUILD_TESTS)
	# Create a set for build all unit tests
	set_property(GLOBAL PROPERTY UNIT_TEST_TARGETS "")
	set_property(GLOBAL PROPERTY UNIT_TEST_SOURCES "")

	get_property(test_targets GLOBAL PROPERTY UNIT_TEST_TARGETS)
	if(test_targets)
		# Create a custom target for all tests
		add_custom_target(all_unit_tests DEPENDS ${test_targets})
	endif()

	get_property(test_source GLOBAL PROPERTY UNIT_TEST_SOURCES)
	if(test_source)
		# Create an executable for the custom target, such that the IDEs can see it as a runnable target
		add_executable(run_all_unit_tests EXCLUDE_FROM_ALL ${test_source})
		# Link the executable with GTest and the TokenValueParser library
		target_link_libraries(run_all_unit_tests PRIVATE GTest::gtest_main GenParseBase GenParseCmdLine)
		set_target_properties(run_all_unit_tests PROPERTIES INSTALLABLE OFF)
	endif()

	set(TESTS_DIRS "")
	foreach(dir ${Shared_DIRS})
		if(IS_DIRECTORY "${dir}/tests")
			list(APPEND TESTS_DIRS "${dir}/tests")
		endif()
	endforeach()

	foreach(dir ${Static_DIRS})
		if(IS_DIRECTORY "${dir}/tests")
			list(APPEND TESTS_DIRS "${dir}/tests")
		endif()
	endforeach()

	#message(STATUS ${TESTS_DIRS})

	# Do not install GTest when packaging Utilities targets
	set(INSTALL_GTEST OFF)
	
	# Add all the tests directories
	foreach(tests_dir ${TESTS_DIRS})
		add_subdirectory("${tests_dir}")
	endforeach()
endif()
# End Create Unit Test Groups *********************************************************************
#**************************************************************************************************