# CMakeList.txt : UltReality::Utilities::Shared

# Create options that are dependent onthis project being top level
option(UTILITIES_SHARED_VERBOSE "Enable verbose messages for UtilitiesShared" ${PROJECT_IS_TOP_LEVEL})

message(STATUS "UTILITIES_SHARED_VERBOSE: ${UTILITIES_SHARED_VERBOSE}")

# Target Creation *********************************************************************************
#**************************************************************************************************

# Create list of all source files to be used in this target
file(GLOB_RECURSE UtilitiesShared_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/*/src/*.cpp")
file(GLOB_RECURSE Shared_DIRS LIST_DIRECTORIES true ${CMAKE_CURRENT_SOURCE_DIR}/*)

set(UtilitiesShared_HEADERS "")
foreach(dir ${Shared_DIRS})
	if(IS_DIRECTORY "${dir}/include")
		list(APPEND UtilitiesShared_HEADERS "${dir}/include")
	endif()
endforeach()
list(APPEND UtilitiesShared_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include")

#message(STATUS ${UtilitiesShared_SOURCE})
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
		TARGETS UtilitiesShared 
		EXPORT UtilitiesShared_Targets 
		ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # Static libraries/import libraries (.lib files for .dll linking) 
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # Shared libraries (.so) 
		RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # .exe or .dll 
		PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # Headers/include directories marked as PUBLIC 
		PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # Headers/include directories marked as PRIVATE
	)

	# Create the targets CMake file which contains the above definitions
	install(
		EXPORT UtilitiesShared_Targets 
		FILE UtilitiesShared_Targets.cmake 
		NAMESPACE UltReality:: 
		DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/UtilitiesShared/cmake
	)

	# Install the actual includes
	install(
		DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ 
		DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	)

	# Generate and install the package version config files
	include(CMakePackageConfigHelpers)
	write_basic_package_version_file(
		"UtilitiesShared_ConfigVersion.cmake" 
		VERSION ${PROJECT_VERSION} 
		COMPATIBILITY SameMajorVersion
	)
	configure_package_config_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/cmake_config/Utilities_Shared_Config.cmake.in" 
		"${CMAKE_CURRENT_BINARY_DIR}/Utilities_Shared_Config.cmake" 
		INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/UtilitiesShared/cmake
	)

	# Install the CMake config files
	install(
		FILES "${CMAKE_CURRENT_BINARY_DIR}/Utilities_Shared_ConfigVersion.cmake" 
		"${CMAKE_CURRENT_BINARY_DIR}/Utilities_Shared_Config.cmake" 
		DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/UtilitiesShared/cmake
	)

	# Define Package install paths
	set(INCLUDEDIR_FOR_PKG_CONFIG "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
	set(LIBDIR_PKG_CONFIG "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

	# Create and install the package config file
	configure_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/cmake_config/Utilities_Shared.pc.in" 
		"${CMAKE_CURRENT_BINARY_DIR}/Utilities_Shared.pc" @ONLY
	)

	# Install the package config file
	install(
		FILES "${CMAKE_CURRENT_BINARY_DIR}/Utilities_Shared.pc" 
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
	)

	# A version that is often used to denote a specific build of the software, including revisions, builds, or other metadata
	set(PACKAGE_VERSION_BUILD "${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}")

	set(PACKAGE_VERSION "${PROJECT_VERSION}-${PACKAGE_VERSION_BUILD}")

	set(CPACK_PACKAGE_DIRECTORY "${CMAKE_SOURCE_DIR}/out/package")

	set(CPACK_PACKAGE_NAME "UtilitiesShared")
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
	set(SHARED_TESTS_DIRS "")

	foreach(dir ${Shared_DIRS})
		if(IS_DIRECTORY "${dir}/tests")
			list(APPEND SHARED_TESTS_DIRS "${dir}/tests")
		endif()
	endforeach()

	#message(STATUS ${TESTS_DIRS})

	# Do not install GTest when packaging Utilities targets
	set(INSTALL_GTEST OFF)
	
	# Add all the tests directories
	foreach(tests_dir ${SHARED_TESTS_DIRS})
		add_subdirectory("${tests_dir}")
	endforeach()
endif()
# End Create Unit Test Groups *********************************************************************
#**************************************************************************************************