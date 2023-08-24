find_package(Sphinx REQUIRED)

#=================================================================== SPHINX ===
# Set sphinx variables
set(SPHINX_HTML_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/html)
set(SPHINX_HTML_FILE ${SPHINX_HTML_OUTPUT_DIR}/index.html)
set(SPHINX_IN ${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in)
set(SPHINX_OUT ${CMAKE_CURRENT_BINARY_DIR}/conf.py)

# Replace variables inside @@ with the current values
configure_file(${SPHINX_IN} ${SPHINX_OUT} @ONLY)

# Explicit list of documentation source files because globbing won't trigger CMake config on re-builds
set(DOCUMENTATION_SOURCE
    ${CMAKE_CURRENT_SOURCE_DIR}/api.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/build.txt
    ${CMAKE_CURRENT_SOURCE_DIR}/changelog.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/contribution.txt
    ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
    ${CMAKE_CURRENT_SOURCE_DIR}/dependencies.txt
    ${CMAKE_CURRENT_SOURCE_DIR}/devops.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/gitlab_setup.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/index.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/md2rst.txt
    ${CMAKE_CURRENT_SOURCE_DIR}/project_brief.txt
    ${CMAKE_CURRENT_SOURCE_DIR}/release_philosophy.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/software.bib
    ${CMAKE_CURRENT_SOURCE_DIR}/targets.txt
    ${CMAKE_CURRENT_SOURCE_DIR}/theory.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/user.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/zreferences.rst
)

# Add sphinx command
add_custom_command(OUTPUT ${SPHINX_HTML_FILE}
                   COMMAND
                       ${CMAKE_COMMAND} -E copy_directory
                           ${CMAKE_CURRENT_SOURCE_DIR}
                           ${CMAKE_CURRENT_BINARY_DIR}
                   COMMAND
                       ${SPHINX_EXECUTABLE} -M html
                       ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   DEPENDS
                       ${DOXYGEN_INDEX_FILE}
                       # Other docs files you want to track should go here (or in some variable)
                       ${DOCUMENTATION_SOURCE}
                   MAIN_DEPENDENCY ${SPHINX_IN} ${SPHINX_OUT}
                   COMMENT "Generating ${PROJECT_NAME} Sphinx documentation")

# Add sphinx target
