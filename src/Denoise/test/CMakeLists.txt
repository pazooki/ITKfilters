include_directories(SYSTEM ${GMOCK_INCLUDE_DIR} ${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Common)
include_directories("..")
# add_library(lpoTest prog_options_test.h)
# target_link_libraries(lpoTest ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(test_denoise.exe test_denoise.cpp)
set_source_files_properties (test_denoise.exe PROPERTIES COMPILE_FLAGS "-Wno-sign-compare")
target_link_libraries(test_denoise.exe ldenoise)
# target_link_libraries(test_denoise.exe ${GMOCK_LIBRARY_MAIN} ${GTEST_LIBRARY_MAIN} ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries(test_denoise.exe gtest_main ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(test_denoise.exe ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_custom_command(
    TARGET test_denoise.exe PRE_BUILD
    COMMAND mkdir -p ${EXECUTABLE_OUTPUT_PATH}/testResults
   )

