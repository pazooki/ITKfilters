include_directories(SYSTEM ${GMOCK_INCLUDE_DIR} ${GTEST_INCLUDE_DIR})
include_directories("..")
file(COPY "../../fixtures" DESTINATION ${EXECUTABLE_OUTPUT_PATH} )
add_executable(test_denoise.exe test_denoise.cpp)
set_source_files_properties (test_denoise.exe PROPERTIES COMPILE_FLAGS "-Wno-sign-compare")
target_link_libraries(test_denoise.exe ldenoise)
target_link_libraries(test_denoise.exe ${GMOCK_LIBRARY_MAIN} ${GTEST_LIBRARY_MAIN} ${CMAKE_THREAD_LIBS_INIT} )


