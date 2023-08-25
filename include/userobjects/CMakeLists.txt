if(xdmf_FOUND)
    target_sources(${PROJECT_NAME} PRIVATE "OutputXdmf.h")
endif()

if(overlap_coupling_FOUND)
    target_sources(${PROJECT_NAME} PRIVATE "OverlapCoupling.h")
endif()

target_sources(${PROJECT_NAME} PRIVATE "ElementIntegrateUserObject.h")
