SET( ITK_LIBRARIES ${ITK_LIBRARIES} PARENT_SCOPE )
include_directories(SYSTEM ${GMOCK_INCLUDE_DIR} ${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Common)
include_directories("..")
# add_library(lpoTest prog_options_test.h)
# target_link_libraries(lpoTest ${Boost_PROGRAM_OPTIONS_LIBRARY})


# test_quadrature
add_executable(test_quadrature.exe test_quadrature.cpp)
set_source_files_properties (test_quadrature.exe PROPERTIES COMPILE_FLAGS "-Wno-sign-compare")
target_link_libraries(test_quadrature.exe gtest_main ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(test_quadrature.exe ${ITK_LIBRARIES})
target_link_libraries(test_quadrature.exe ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_custom_command(
    TARGET test_quadrature.exe PRE_BUILD
    COMMAND mkdir -p ${EXECUTABLE_OUTPUT_PATH}/testQuadrature
   )
