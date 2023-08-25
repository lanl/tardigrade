if(xdmf_FOUND)
    target_sources(${PROJECT_NAME} PRIVATE "OutputXdmf.C")
endif()

if(overlap_coupling_FOUND)
    target_sources(${PROJECT_NAME} PRIVATE "OverlapCoupling.C")
endif()

target_sources(${PROJECT_NAME} PRIVATE "ElementIntegrateUserObject.C")
