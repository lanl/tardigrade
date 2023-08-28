set(${PROJECT_NAME}_userobjects_TESTS "overlapping_columns")

set(${PROJECT_NAME}_overlap_coupling_TESTS "overlapping_columns")

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

foreach(filename ${${PROJECT_NAME}_userobjects_TESTS})
    if(filename IN_LIST EXCLUDE_FILENAMES)
    else()
        set(${PROJECT_NAME}_regression_TESTS ${${PROJECT_NAME}_regression_TESTS} ${filename})
    endif()
endforeach(filename)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    foreach(test ${${PROJECT_NAME}_regression_TESTS})
        add_subdirectory(${test})
    endforeach(test)
endif()
