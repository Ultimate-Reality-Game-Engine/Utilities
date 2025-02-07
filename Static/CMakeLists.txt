# CMakeList.txt : UltReality::Utilities::Static

# Target Creation *********************************************************************************
#**************************************************************************************************

file(GLOB_RECURSE UtilitiesStatic_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/*/src/*.cpp")
file(GLOB_RECURSE Static_DIRS LIST_DIRECTORIES true ${CMAKE_CURRENT_SOURCE_DIR}/*)

set(UtilitiesStatic_HEADERS "")
foreach(dir ${Static_DIRS})
	if(IS_DIRECTORY "${dir}/include")
		list(APPEND UtilitiesStatic_HEADERS "${dir}/include")
	endif()

	if(IS_DIRECTORY "${dir}/inl")
		list(APPEND UtilitiesStatic_HEADERS "${dir}/inl")
	endif()
endforeach()

#message(STATUS ${UtilitiesStatic_SOURCE})
#message(STATUS ${UtilitiesStatic_HEADERS})

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
		#target_link_libraries(UtilitiesStatic )

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
	# Install the targets
	install(
		TARGETS UtilitiesStatic
		EXPORT UtilitiesStatic_Targets 
		ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # Static libraries/import libraries (.lib files for .dll linking) 
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # Shared libraries (.so) 
		RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # .exe or .dll 
		PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # Headers/include directories marked as PUBLIC 
		PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # Headers/include directories marked as PRIVATE
	)

	# Create the targets CMake file which contains the above definitions
	install(
		EXPORT UtilitiesStatic_Targets 
		FILE UtilitiesStatic_Targets.cmake 
		NAMESPACE UltReality:: 
		DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/UtilitiesStatic/cmake
	)

	# Install the actual includes
	install(
		DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ 
		DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	)

	# Generate and install the package version config files
	include(CMakePackageConfigHelpers)
	write_basic_package_version_file(
		"UtilitiesStatic_ConfigVersion.cmake" 
		VERSION ${PROJECT_VERSION} 
		COMPATIBILITY SameMajorVersion
	)
	configure_package_config_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/cmake_config/Utilities_Static_Config.cmake.in" 
		"${CMAKE_CURRENT_BINARY_DIR}/Utilities_Static_Config.cmake" 
		INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/UtilitiesStatic/cmake
	)

	# Install the CMake config files
	install(
		FILES "${CMAKE_CURRENT_BINARY_DIR}/Utilities_Static_ConfigVersion.cmake" 
		"${CMAKE_CURRENT_BINARY_DIR}/Utilities_Static_Config.cmake" 
		DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/UtilitiesStatic/cmake
	)

	# Define Package install paths
	set(INCLUDEDIR_FOR_PKG_CONFIG "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
	set(LIBDIR_PKG_CONFIG "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

	# Create and install the package config file
	configure_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/cmake_config/Utilities_Static.pc.in" 
		"${CMAKE_CURRENT_BINARY_DIR}/Utilities_Static.pc" @ONLY
	)

	# Install the package config file
	install(
		FILES "${CMAKE_CURRENT_BINARY_DIR}/Utilities_Static.pc" 
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
	)

	# A version that is often used to denote a specific build of the software, including revisions, builds, or other metadata
	set(PACKAGE_VERSION_BUILD "${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}")

	set(PACKAGE_VERSION "${PROJECT_VERSION}-${PACKAGE_VERSION_BUILD}")

	set(CPACK_PACKAGE_DIRECTORY "${CMAKE_SOURCE_DIR}/out/package")

	set(CPACK_PACKAGE_NAME "UtilitiesStatic")
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
	set(STATIC_TESTS_DIRS "")

	foreach(dir ${Static_DIRS})
		if(IS_DIRECTORY "${dir}/tests")
			list(APPEND STATIC_TESTS_DIRS "${dir}/tests")
		endif()
	endforeach()

	#message(STATUS ${TESTS_DIRS})

	# Do not install GTest when packaging Utilities targets
	set(INSTALL_GTEST OFF)
	
	# Add all the tests directories
	foreach(tests_dir ${STATIC_TESTS_DIRS})
		add_subdirectory("${tests_dir}")
	endforeach()
endif()
# End Create Unit Test Groups *********************************************************************
#**************************************************************************************************