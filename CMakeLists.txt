#===================================================================================================== PROJECT SETUP ===
cmake_minimum_required(VERSION 3.18)
project(tardigrade)
set(PROJECT_EXECUTABLE_NAME "tardigrade-opt")

# Set common project paths relative to project root directory
set(CPP_SRC_PATH "src")
set(CPP_INC_PATH "include")
set(CPP_TEST_PATH "test")
set(CPP_REGRESSION_TEST_PATH "test/tests")
set(CMAKE_SRC_PATH "src/cmake")
set(DOXYGEN_SRC_PATH "docs/doxygen")
set(SPHINX_SRC_PATH "docs/sphinx")

# Add the cmake folder to locate project CMake module(s)
set(XDMF_INSTALL_DIR "/path/to/Xdmf-config.cmake" CACHE PATH "The location of the Xdmf-config.cmake or XdmfConfig.cmake file")
set(VORO_INSTALL_DIR "/path/to/voro-config.cmake" CACHE PATH "The location of the voro-config.cmake or VOROConfig.cmake file")

set(ADDITIONAL_INSTALL_DIRECTORIES ${XDMF_INSTALL_DIR} ${VORO_INSTALL_DIR})
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${ADDITIONAL_INSTALL_DIRECTORIES})

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}" ${CMAKE_MODULE_PATH})

# Set build type checks
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_lower)
set(upstream_required "")
set(not_conda_test "true")
if(cmake_build_type_lower STREQUAL "release")
    set(upstream_required "REQUIRED")
elseif(cmake_build_type_lower STREQUAL "conda-test")
    set(upstream_required "REQUIRED")
    set(not_conda_test "false")
    # Find the installed project libraries
    find_package(${PROJECT_NAME} ${upstream_required} CONFIG)
    # Find the installed project umat
endif()

# Get version number from Git
set(VERSION_UPDATE_FROM_GIT True)
if(${not_conda_test} STREQUAL "true")
    # TODO: On osx-arm64 some c++ projects struggle to find $PREFIX/bin/python during conda-builds and others don't.
    # Trace the source and remove the hint when CMake configures with the correct python version during conda-build.
    if (DEFINED ENV{PREFIX})
        set(Python_ROOT_DIR "$ENV{PREFIX}/bin")
    endif()
    find_package(Python COMPONENTS Interpreter REQUIRED)
    execute_process(COMMAND ${Python_EXECUTABLE} -m setuptools_scm
                            OUTPUT_VARIABLE ${PROJECT_NAME}_VERSION_STRING_FULL
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX MATCH "[0-9]+\.[0-9]+\.[0-9]+" ${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_VERSION_STRING_FULL})
    project(${PROJECT_NAME} VERSION ${${PROJECT_NAME}_VERSION})
endif()

# Add installation directory variables
include(GNUInstallDirs)

# Make the code position independent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set the c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt")
endif()

# Enable CTest
enable_testing()

#================================================================================================= FIND DEPENDENCIES ===
# Find eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
if(EIGEN3_FOUND)
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
endif()

# Find Xdmf
find_package(xdmf NAMES "Xdmf" "xdmf")
include_directories(${XDMF_INCLUDE_DIRS})
link_directories(${XDMF_CMAKE_INSTALL_FULL_LIBDIR})
link_directories("${XDMF_CMAKE_INSTALL_FULL_LIBDIR}/../lib") #XDMF seems to be installing files into lib and lib64
if(xdmf_FOUND)
    message(STATUS "Found Xdmf: ${XDMF_INCLUDE_DIRS}")

    # voro++ is only used if the overlap coupling module is used which can only be built
    # if xdmf is found
    find_package(VORO REQUIRED)
    include_directories(${VORO_INCLUDE_DIR})
    if(VORO_FOUND)
        message(STATUS "Found voro++: ${VORO_INCLUDE_DIR}")
    endif()

else()
    message(WARNING "Xdmf not found. Building without Xdmf components.")
endif()

# Find related, but external, projects in installed environment
set(upstream_packages "tardigrade_error_tools" "tardigrade_vector_tools" "tardigrade_micromorphic_element")
# Add overlap coupling if required. At this time, using overlap_coupling requires
# that xdmf is installed. This requirement may change in time
set(optional_upstream_packages "overlap_coupling")
include(FetchContent)
foreach(package ${upstream_packages})
    string(TOUPPER "${package}" package_upper)
    find_package(${package} ${upstream_required} CONFIG)
    if(${package}_FOUND)
        message(STATUS "Found ${package}: ${${package}_DIR}")
    else()
        # Find related, but external, projects using FetchContent and building locally
        message(WARNING "Did not find an installed ${package} package. Attempting local build with FetchContent.")
        # TODO: Decide how to handle main/dev build differences. bash env variables for GIT TAG?
        FetchContent_Declare(
            ${package}
            GIT_REPOSITORY ssh://git@re-git.lanl.gov:10022/aea/material-models/${package}.git
            GIT_TAG        origin/dev
        )
        FetchContent_MakeAvailable(${package})
    endif()
endforeach(package)

foreach(package ${optional_upstream_packages})
    string(TOUPPER "${package}" package_upper)
    find_package(${package} CONFIG)
    if(${package}_FOUND)
        message(STATUS "Found ${package}: ${${package}_DIR}")
    else()
        # Find related, but external, projects using FetchContent and building locally
        message(WARNING "Did not find an installed ${package} package. This package is optional.")
    endif()
endforeach(package)

if(DEFINED ENV{MOOSE_DIR})
    message(STATUS "MOOSE_DIR defined: $ENV{MOOSE_DIR}")
elseif(EXISTS "${PROJECT_BINARY_DIR}/deps_/moose-src") #Prevents pulling the repo if we've already done that
    message(STATUS "MOOSE_DIR not defined but a previously existing Git fetch discovered")
    set(ENV{MOOSE_DIR} "${PROJECT_BINARY_DIR}/deps_/moose-src")
    message(STATUS "  Setting MOOSE_DIR to: $ENV{MOOSE_DIR}")
else()
    message(STATUS "No Moose directory found. Will attempt to clone.")

    set(MOOSE_REPO_URL "https://github.com/idaholab/moose.git")

    message(STATUS "  Fetching from ${MOOSE_REPO_URL}")
    set(FETCH_CONTENT_QUIET FALSE)
    FetchContent_Declare(
        moose
        GIT_REPOSITORY ${MOOSE_REPO_URL}
        GIT_PROGRESS TRUE
        GIT_TAG origin/master
    )

    FetchContent_MakeAvailable(moose)
    set(ENV{MOOSE_DIR} ${moose_SOURCE_DIR})
    message("  Setting MOOSE_DIR to: $ENV{MOOSE_DIR}")
endif()



## Fetch additional header files
#FetchContent_Declare(
#    quickhull
#    GIT_REPOSITORY https://github.com/akuukka/quickhull.git
#    GIT_TAG        origin/master
#)
#FetchContent_MakeAvailable(quickhull)
#message("quickhull found: ${quickhull_SOURCE_DIR}")
#include_directories(${quickhull_SOURCE_DIR})

##=============================================================================================== ADD PROJECT TARGETS ===
# MUST COME AFTER DEPENDENCY LOCATING
# Add project source directories

# Configure MOOSE for installation
if(cmake_build_type_lower STREQUAL release)
    set(MOOSE_CONFIGURE_COMMAND "$ENV{MOOSE_DIR}/configure")
    set(MOOSE_CONFIGURE_OPTIONS "--prefix=${MOOSE_INSTALL_PREFIX}")
    execute_process(COMMAND ${MOOSE_CONFIGURE_COMMAND} ${MOOSE_CONFIGURE_OPTIONS}
                    WORKING_DIRECTORY $ENV{MOOSE_DIR}
                    RESULT_VARIABLE MOOSE_CONFIGURE_RESULT
                    OUTPUT_VARIABLE MOOSE_CONFIGURE_OUTPUT)
    message("MOOSE_CONFIGURE_RESULT:\n${MOOSE_CONFIGURE_RESULT}")
    message("MOOSE_CONFIGURE_OUTPUT:\n${MOOSE_CONFIGURE_OUTPUT}")
endif()

add_custom_target(${PROJECT_NAME} ALL
                  COMMAND MOOSE_DIR=$ENV{MOOSE_DIR} make -j4 -f "${PROJECT_NAME}_Makefile"
                  WORKING_DIRECTORY ${PROJECT_BINARY_DIRECTORY})

include_directories(${CPP_INC_PATH})
add_subdirectory(${CPP_SRC_PATH})
add_subdirectory(${CPP_INC_PATH})
add_subdirectory(${CPP_TEST_PATH})

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)
    add_subdirectory(${CPP_REGRESSION_TEST_PATH})
    target_sources(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/testroot")
    target_sources(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/run_tests")
    add_test(NAME "${PROJECT_NAME}_Regression_Tests"
             COMMAND run_tests
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
    set_tests_properties("${PROJECT_NAME}_Regression_Tests" PROPERTIES
                         ENVIRONMENT MOOSE_DIR=$ENV{MOOSE_DIR} )
    # Add docs
    if(${not_conda_test} STREQUAL "true")
        add_subdirectory(${DOXYGEN_SRC_PATH})
        add_subdirectory(${SPHINX_SRC_PATH})
    endif()
endif()


# Configure the tardigrade makefile
configure_file("${PROJECT_SOURCE_DIR}/Makefile.in" "${PROJECT_BINARY_DIR}/${PROJECT_NAME}_Makefile")

# Copy all of the projects source files to the build directory
get_target_property(${PROJECT_NAME}_SOURCES ${PROJECT_NAME} SOURCES)

foreach(source_file ${${PROJECT_NAME}_SOURCES})
    file(RELATIVE_PATH rel ${PROJECT_SOURCE_DIR} ${source_file})
    configure_file("${source_file}" "${PROJECT_BINARY_DIR}/${rel}" COPYONLY)
endforeach(source_file)

#==================================================================================== SETUP INSTALLATION CMAKE FILES ===
foreach(package ${PROJECT_NAME})
    include(CMakePackageConfigHelpers)
    write_basic_package_version_file("${package}ConfigVersion.cmake"
                                     VERSION ${PROJECT_VERSION}
                                     COMPATIBILITY SameMajorVersion)
    configure_package_config_file(
      "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}/Config.cmake.in"
      "${PROJECT_BINARY_DIR}/${package}Config.cmake"
      INSTALL_DESTINATION
      ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)

    # CMake won't build the targets for local builds of upstream projects
    if(cmake_build_type_lower STREQUAL release)

#        set(${PROJECT_NAME}_INSTALL_COMMAND "execute_process(COMMAND MOOSE_DIR=$ENV{MOOSE_DIR} make -f ${PROJECT_NAME}_Makefile install WORKING_DIRECTORY ${PROJECT_BINARY_DIR} RESULT_VARIABLE tmp_res OUTPUT_VARIABLE tmp_out)")
        set(${PROJECT_NAME}_INSTALL_COMMAND "execute_process(COMMAND /bin/bash -c \"MOOSE_DIR=$ENV{MOOSE_DIR} make -f ${PROJECT_NAME}_Makefile install\" WORKING_DIRECTORY ${PROJECT_BINARY_DIR} RESULT_VARIABLE tmp_res OUTPUT_VARIABLE tmp_out)")
        message("${PROJECT_NAME}_INSTALL_COMMAND: ${${PROJECT_NAME}_INSTALL_COMMAND}")
    
        install(CODE "${${PROJECT_NAME}_INSTALL_COMMAND}")
        install(CODE [=[message("MOOSE INSTALL ERRORS:\n${tmp_res}")]=])
        install(CODE [=[message("MOOSE INSTALL OUTPUT:\n${tmp_out}")]=])

    endif()

    install(FILES "${PROJECT_BINARY_DIR}/${package}Config.cmake"
                  "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
            DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)


endforeach(package)
