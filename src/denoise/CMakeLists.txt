SET( OpenCV_LIBS ${OpenCV_LIBS} PARENT_SCOPE )
SET( ITK_LIBRARIES ${ITK_LIBRARIES} PARENT_SCOPE )
SET( Boost_LIBRARIES ${Boost_LIBRARIES} PARENT_SCOPE )

add_library(ldenoise STATIC ldenoise.cpp)
target_link_libraries(SAXSsim ${Boost_LIBRARIES} ${OpenCV_LIBS})
target_link_libraries(SAXSsim ${ITK_LIBRARIES} )

message(status "** OpenCV Libraries: ${OpenCV_LIBS}")
message(status "** ITK Libraries: ${ITK_LIBRARIES}")
message(status "** Boost Libraries: ${Boost_LIBRARIES}")

SET(ENABLE_GOOGLE_TEST ${ENABLE_GOOGLE_TEST} PARENT_SCOPE)
IF(ENABLE_GOOGLE_TEST)
    add_subdirectory(test)
ENDIF()
