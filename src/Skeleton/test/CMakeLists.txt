SET( ITK_LIBRARIES ${ITK_LIBRARIES} PARENT_SCOPE )
SET( Boost_LIBRARIES ${Boost_LIBRARIES} PARENT_SCOPE )
include_directories(SYSTEM ${GMOCK_INCLUDE_DIR} ${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Common)
include_directories("..")
add_executable(test_skeleton.exe test_skeleton.cpp ../Skeleton.h)
# set_source_files_properties (test_skeleton.exe PROPERTIES COMPILE_FLAGS "-Wno-sign-compare")
# target_link_libraries(test_skeleton.exe ${GMOCK_LIBRARY_MAIN} ${GTEST_LIBRARY_MAIN} ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries(test_skeleton.exe gtest_main ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(test_skeleton.exe ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_link_libraries(test_skeleton.exe ${ITK_LIBRARIES})

add_custom_command(
    TARGET test_skeleton.exe PRE_BUILD
    COMMAND mkdir -p ${EXECUTABLE_OUTPUT_PATH}/testResults
   )

