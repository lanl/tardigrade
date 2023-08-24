set(${PROJECT_NAME}_regression_tests_DIRECTORIES "kernels" "userobjects")

set(${PROJECT_NAME}_kernels_TESTS "const_mms" "linear_mms" "rotate_stretch" "stretch_y" "const_mms_overlap" "linear_mms_overlap" "simple_diffusion" "stretch_y_plastic")

set(${PROJECT_NAME}_userobjects_TESTS "overlapping_columns")

set(${PROJECT_NAME}_overlap_coupling_TESTS "kernels/const_mms_overlap" "kernels/linear_mms_overlap" "userobjects/overlapping_columns")

set(EXCLUDE_FILENAMES "")
foreach(package ${optional_upstream_packages})
    if(${${package}_FOUND})
    else()
        foreach(filename "${${PROJECT_NAME}_${package}_TESTS}")
            set(EXCLUDE_FILENAMES ${EXCLUDE_FILENAMES} ${filename})
        endforeach(filename)
    endif()
endforeach(package)

set(${PROJECT_NAME}_regression_TESTS "")

foreach(directory ${${PROJECT_NAME}_regression_tests_DIRECTORIES})
    foreach(filename ${${PROJECT_NAME}_${directory}_TESTS})
        set(root_path "${directory}/${filename}")
        if(${root_path} IN_LIST EXCLUDE_FILENAMES)
        else()
            set(${PROJECT_NAME}_regression_TESTS ${${PROJECT_NAME}_regression_TESTS} ${root_path})
        endif()
    endforeach(filename)
endforeach(directory)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    foreach(test ${${PROJECT_NAME}_regression_TESTS})
        add_subdirectory(${test})
    endforeach(test)
endif()
