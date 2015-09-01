SET( ITK_LIBRARIES ${ITK_LIBRARIES} PARENT_SCOPE )
SET( Boost_LIBRARIES ${Boost_LIBRARIES} PARENT_SCOPE )

# add_library(lskeleton Skeleton.cpp)
# target_link_libraries(lskeleton ${Boost_LIBRARIES})
# target_link_libraries(lskeleton ${ITK_LIBRARIES} )

# add_executable(skeleton.exe Skeleton.cpp)
# target_link_libraries(skeleton.exe ${ITK_LIBRARIES} )

message(status "** ITK Libraries: ${ITK_LIBRARIES}")
message(status "** Boost Libraries: ${Boost_LIBRARIES}")

SET(ENABLE_GOOGLE_TEST ${ENABLE_GOOGLE_TEST} PARENT_SCOPE)
IF(ENABLE_GOOGLE_TEST)
    add_subdirectory(test)
ENDIF()
