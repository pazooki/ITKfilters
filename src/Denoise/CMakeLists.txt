SET( ITK_LIBRARIES ${ITK_LIBRARIES} PARENT_SCOPE )
SET( Boost_LIBRARIES ${Boost_LIBRARIES} PARENT_SCOPE )

add_library(ldenoise STATIC Denoise.cpp)
target_link_libraries(ldenoise ${Boost_LIBRARIES})
target_link_libraries(ldenoise ${ITK_LIBRARIES} )
add_library(lFlatStructuringElementWithImageBridge STATIC itkFlatStructuringElementWithImageBridge.h itkFlatStructuringElementWithImageBridge.hxx)
target_link_libraries(lFlatStructuringElementWithImageBridge ${ITK_LIBRARIES} )

message(status "** ITK Libraries: ${ITK_LIBRARIES}")
message(status "** Boost Libraries: ${Boost_LIBRARIES}")

SET(ENABLE_GOOGLE_TEST ${ENABLE_GOOGLE_TEST} PARENT_SCOPE)
IF(ENABLE_GOOGLE_TEST)
    add_subdirectory(test)
    add_subdirectory(tem_saxs_paper)
ENDIF()
