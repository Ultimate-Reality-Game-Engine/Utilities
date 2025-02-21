# CMakeList.txt : UltReality::Utilities::Shared

project(UTILITIES_SHARED VERSION 1.0.0)

# Create options that are dependent onthis project being top level
option(UTILITIES_SHARED_VERBOSE "Enable verbose messages for UtilitiesShared" ${PROJECT_IS_TOP_LEVEL})

message(STATUS "UTILITIES_SHARED_VERBOSE: ${UTILITIES_SHARED_VERBOSE}")

# Target Creation *********************************************************************************
#**************************************************************************************************

option(UTILITIES_SHARED_DEBUG "Enable CMake related debug messages" OFF)

# Create list of all source files to be used in this target
file(GLOB_RECURSE UtilitiesShared_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/*/src/*.cpp")
file(GLOB Shared_DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*)

#message(STATUS "Shared_DIRS: ${Shared_DIRS}")

if(NOT DEFINED UTILITIES_SHARED_BUILD)
	set(UTILITIES_SHARED_BUILD ON)
endif()

if(UTILITIES_SHARED_BUILD)
	# Create the UtilitiesShared target
	if(UtilitiesShared_SOURCE)
		add_library(UtilitiesShared SHARED ${UtilitiesShared_SOURCE})

		target_compile_definitions(UtilitiesShared PRIVATE LIBRARY_EXPORTS)

		# Set up include directories for the shared library target
		target_include_directories(UtilitiesShared PUBLIC 
			$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>  # Add Shared/include explicitly
			$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/UtilitiesShared> # Install it at the top-level include
		)

		foreach(dir ${Shared_DIRS})
			if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/include")
				if(UTILITIES_SHARED_DEBUG)	
					message(STATUS "Adding include directory: ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/include")
				endif()
				target_include_directories(UtilitiesShared PUBLIC 
					$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${dir}/include> 
 					$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/UtilitiesShared/${dir}> # This is used when the library is installed
				)
			endif()

			if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/inl")
				if(UTILITIES_SHARED_DEBUG)
					message(STATUS "Adding inl directory: ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/inl")
				endif()
				target_include_directories(UtilitiesShared PUBLIC 
					$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${dir}/inl> 
 					$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/UtilitiesShared/${dir}> # This is used when the library is installed
				)
			endif()
		endforeach()

		# Link libraries to the target
		#target_link_libraries()

		set_target_properties(UtilitiesShared PROPERTIES 
			VERSION ${PROJECT_VERSION} 
			SOVERSION ${PROJECT_VERSION_MAJOR}
		)

# End Target Creation *****************************************************************************
#**************************************************************************************************

# Installation and Packing Configuration **********************************************************
#**************************************************************************************************

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
			NAMESPACE UltReality::Utilities 
			DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/UtilitiesShared
		)

		# Install the actual includes
		install(
			DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include"
			DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/UtilitiesShared"
		)

		foreach(dir ${Shared_DIRS})
			if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/include")
				install(
					DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/include/"
					DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/UtilitiesShared/${dir}"
				)
			endif()

			if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/inl")
				install(
					DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/inl/"
					DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/UtilitiesShared/${dir}"
				)
			endif()
		endforeach()

		# Generate and install the package version config files
		include(CMakePackageConfigHelpers)
		write_basic_package_version_file(
			"Utilities_Shared_ConfigVersion.cmake" 
			VERSION ${PROJECT_VERSION} 
			COMPATIBILITY SameMajorVersion
		)
		configure_package_config_file(
			"${CMAKE_CURRENT_SOURCE_DIR}/cmake_config/Utilities_Shared_Config.cmake.in" 
			"${CMAKE_CURRENT_BINARY_DIR}/Utilities_Shared_Config.cmake" 
			INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/UtilitiesShared
		)

		# Install the CMake config files
		install(
			FILES "${CMAKE_CURRENT_BINARY_DIR}/Utilities_Shared_ConfigVersion.cmake" 
			"${CMAKE_CURRENT_BINARY_DIR}/Utilities_Shared_Config.cmake" 
			DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/UtilitiesShared
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

		if(UTILITIES_SHARED_VERBOSE)
			message(STATUS "PACKAGE_VERSION is: ${PACKAGE_VERSION}")
			message(STATUS "PACKAGE_FILE_NAME is: ${CPACK_PACKAGE_FILE_NAME}")
		endif()

# End Installation and Packing Configuration ******************************************************
#**************************************************************************************************

# Create Unit Test Groups *************************************************************************
#**************************************************************************************************
		if (UTILITIES_BUILD_TESTS)
			if(UTILITIES_SHARED_DEBUG)
				message(STATUS "Building test suit for Utilities Shared")
			endif()

			set(SHARED_TEST_DIRS "")

			foreach(dir ${Static_DIRS})
				if(IS_DIRECTORY "${dir}/tests")
					list(APPEND SHARED_TEST_DIRS "${dir}/tests")
				endif()
				if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/tests")
						if(UTILITIES_SHARED_DEBUG)	
							message(STATUS "Adding test directory: ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/tests")
						endif()
						list(APPEND SHARED_TEST_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/tests")
					endif()
			endforeach()

			# Do not install GTest when packaging Utilities targets
			set(INSTALL_GTEST OFF)
	
			# Add all the tests directories
			foreach(tests_dir ${SHARED_TEST_DIRS})
				if(UTILITIES_SHARED_DEBUG)
					message(STATUS "Adding Sub-Directory: ${tests_dir}")
				endif()
				add_subdirectory("${tests_dir}")
			endforeach()
		endif()
	endif()
endif()
# End Create Unit Test Groups *********************************************************************
#**************************************************************************************************