set(tardigrade_test_SRC_FILES "")
set(tardigrade_test_INC_FILES "")

set(tardigrade_test_FILE_DIRECTORIES "base")

set(test_base_FILENAMES "tardigradeTestApp")

set(EXCLUDE_FILENAMES "")

foreach(source_directory ${tardigrade_test_FILE_DIRECTORIES})
    foreach(filename ${test_${source_directory}_FILENAMES})
        if ("${filename}" IN_LIST EXCLUDE_FILENAMES) 
        else()
            set(root_inc_dir "${CPP_INC_PATH}")
            set(root_src_dir "${CPP_SRC_PATH}")
            set(current_inc_path "${root_inc_dir}/${source_directory}/${filename}.h")
            set(current_src_path "${root_src_dir}/${source_directory}/${filename}.C")
            set(tardigrade_test_INC_FILES ${tardigrade_INC_FILES} ${current_inc_path})
            set(tardigrade_test_SRC_FILES ${tardigrade_SRC_FILES} ${current_src_path})
        endif()
    endforeach(filename)
endforeach(source_directory)

foreach(src_file ${tardigrade_test_SRC_FILES})
    configure_file(${src_file} "${PROJECT_BINARY_DIR}/${PROJECT_NAME}/${CPP_TEST_PATH}/${src_file}" COPYONLY)
endforeach(src_file)

foreach(inc_file ${tardigrade_test_INC_FILES})
    configure_file(${inc_file} "${PROJECT_BINARY_DIR}/${PROJECT_NAME}/${CPP_TEST_PATH}/${inc_file}" COPYONLY)
endforeach(inc_file)

