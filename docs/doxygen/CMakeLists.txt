find_package(Doxygen REQUIRED)

# Find all the public headers
file(GLOB_RECURSE PROJECT_PUBLIC_HEADERS ${PROJECT_SOURCE_DIR}/${CPP_SRC_PATH}/*.h)
file(GLOB_RECURSE PROJECT_SOURCE_FILES ${PROJECT_SOURCE_DIR}/${CPP_SRC_PATH}/*.cpp)

#=================================================================== DOXYGEN ===
# Set doxygen variables
set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/${CPP_SRC_PATH})  # required by Doxyfile.in
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})  # required by Doxyfile.in
set(DOXYGEN_INDEX_FILE ${CMAKE_CURRENT_BINARY_DIR}/xml/index.xml)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

# Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

# Create doxygen output and command
add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS
                       ${PROJECT_PUBLIC_HEADERS}
                       ${PROJECT_SOURCE_FILES}
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Generating ${PROJECT_NAME} Doxygen API"
                   VERBATIM)

# Add doxygen target
add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})
