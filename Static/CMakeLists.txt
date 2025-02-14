# CMakeList.txt : UltReality::Utilities::Static

project(UTILITIES_STATIC VERSION 1.0.0)

# Create options that are dependent onthis project being top level
option(UTILITIES_STATIC_VERBOSE "Enable verbose messages for UtilitiesStatic" ${PROJECT_IS_TOP_LEVEL})

message(STATUS "UTILITIES_STATIC_VERBOSE: ${UTILITIES_STATIC_VERBOSE}")

# Target Creation *********************************************************************************
#**************************************************************************************************

option(UTILITIES_STATIC_DEBUG "Enable CMake related debug messages" ON)

file(GLOB_RECURSE UtilitiesStatic_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/*/src/*.cpp")
file(GLOB Static_DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*)

#message(STATUS "Static_DIRS: ${Static_DIRS}")
#message(STATUS ${UtilitiesStatic_SOURCE})

if(NOT DEFINED UTILITIES_STATIC_BUILD)
	set(UTILITIES_STATIC_BUILD ON)
endif()

if(UTILITIES_STATIC_BUILD)
	# Create the UtilitiesStatic target
	if(UtilitiesStatic_SOURCE)
		add_library(UtilitiesStatic STATIC ${UtilitiesStatic_SOURCE})

		# Set up include directories for the shared library target
		foreach(dir ${Static_DIRS})
			if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/include")
				if(UTILITIES_STATIC_DEBUG)	
					message(STATUS "Adding include directory: ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/include")
				endif()
				target_include_directories(UtilitiesStatic PUBLIC 
					$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${dir}/include> 
 					$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/UtilitiesStatic/${dir}> # This is used when the library is installed
				)
			endif()

			if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/inl")
				if(UTILITIES_STATIC_DEBUG)	
					message(STATUS "Adding inl directory: ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/inl")
				endif()
				target_include_directories(UtilitiesStatic PUBLIC 
					$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${dir}/inl> 
 					$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/UtilitiesStatic/${dir}> # This is used when the library is installed
				)
			endif()
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
	DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/UtilitiesStatic
)

# Install the actual includes
foreach(dir ${Static_DIRS})
    if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/include")
        install(
            DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/include/"
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/UtilitiesStatic/${dir}"
        )
    endif()

	if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/inl")
        install(
            DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/inl/"
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/UtilitiesStatic/${dir}"
        )
    endif()
endforeach()

# Generate and install the package version config files
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	"Utilities_Static_ConfigVersion.cmake" 
	VERSION ${PROJECT_VERSION} 
	COMPATIBILITY SameMajorVersion
)
configure_package_config_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake_config/Utilities_Static_Config.cmake.in" 
	"${CMAKE_CURRENT_BINARY_DIR}/Utilities_Static_Config.cmake" 
	INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/UtilitiesStatic
)

# Install the CMake config files
install(
	FILES "${CMAKE_CURRENT_BINARY_DIR}/Utilities_Static_ConfigVersion.cmake" 
	"${CMAKE_CURRENT_BINARY_DIR}/Utilities_Static_Config.cmake" 
	DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/UtilitiesStatic
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

if(UTILITIES_STATIC_VERBOSE)
	message(STATUS "PACKAGE_VERSION is: ${PACKAGE_VERSION}")
	message(STATUS "PACKAGE_FILE_NAME is: ${CPACK_PACKAGE_FILE_NAME}")
endif()

# End Installation and Packing Configuration ******************************************************
#**************************************************************************************************

# Create Unit Test Groups *************************************************************************
#**************************************************************************************************
if (UTILITIES_BUILD_TESTS)
	if(UTILITIES_STATIC_DEBUG)
		message(STATUS "Building test suit for Utilities Static")
	endif()

	set(STATIC_TEST_DIRS "")

	foreach(dir ${Static_DIRS})
		if(IS_DIRECTORY "${dir}/tests")
			list(APPEND STATIC_TEST_DIRS "${dir}/tests")
		endif()
		if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/tests")
				if(UTILITIES_STATIC_DEBUG)	
					message(STATUS "Adding test directory: ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/tests")
				endif()
				list(APPEND STATIC_TEST_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/tests")
			endif()
	endforeach()

	# Do not install GTest when packaging Utilities targets
	set(INSTALL_GTEST OFF)
	
	# Add all the tests directories
	foreach(tests_dir ${STATIC_TEST_DIRS})
		if(UTILITIES_STATIC_DEBUG)
			message(STATUS "Adding Sub-Directory: ${tests_dir}")
		endif()
		add_subdirectory("${tests_dir}")
	endforeach()
endif()
# End Create Unit Test Groups *********************************************************************
#**************************************************************************************************