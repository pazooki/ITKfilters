SET( ITK_LIBRARIES ${ITK_LIBRARIES} PARENT_SCOPE )
include_directories(SYSTEM ${GMOCK_INCLUDE_DIR} ${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Common)
include_directories("..")
# add_library(lpoTest prog_options_test.h)
# target_link_libraries(lpoTest ${Boost_PROGRAM_OPTIONS_LIBRARY})
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
find_package(Qt5 COMPONENTS Widgets)

# # test_quadrature
# add_executable(test_quadrature.exe test_quadrature.cpp)
# set_source_files_properties (test_quadrature.exe PROPERTIES COMPILE_FLAGS "-Wno-sign-compare")
# target_link_libraries(test_quadrature.exe gtest_main ${CMAKE_THREAD_LIBS_INIT})
# target_link_libraries(test_quadrature.exe ${ITK_LIBRARIES})
# target_link_libraries(test_quadrature.exe ${Boost_PROGRAM_OPTIONS_LIBRARY})
# target_link_libraries(test_quadrature.exe ${VTK_LIBRARIES})
# target_link_libraries(test_quadrature.exe Qt5::Widgets)
#
# test_riesz
add_executable(test_riesz.exe test_riesz.cpp)
set_source_files_properties (test_riesz.exe PROPERTIES COMPILE_FLAGS "-Wno-sign-compare")
target_link_libraries(test_riesz.exe gtest_main ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(test_riesz.exe ${ITK_LIBRARIES})
target_link_libraries(test_riesz.exe ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_link_libraries(test_riesz.exe ${VTK_LIBRARIES})
target_link_libraries(test_riesz.exe Qt5::Widgets)

# test_wavelet
# add_executable(test_wavelet.exe test_wavelet.cpp /home/phc/Software/ITK/ITK-development/Modules/Filtering/Wavelet/src/itkWaveletGenerator.cxx)
# set_source_files_properties (test_wavelet.exe PROPERTIES COMPILE_FLAGS "-Wno-sign-compare")
# target_link_libraries(test_wavelet.exe gtest_main ${CMAKE_THREAD_LIBS_INIT})
# target_link_libraries(test_wavelet.exe ${ITK_LIBRARIES})
# target_link_libraries(test_wavelet.exe ${Boost_PROGRAM_OPTIONS_LIBRARY})
# target_link_libraries(test_wavelet.exe ${VTK_LIBRARIES})
# target_link_libraries(test_wavelet.exe Qt5::Widgets)
# message( ${ITK_LIBRARIES})
# add_custom_command(
#     TARGET test_quadrature.exe PRE_BUILD
#     COMMAND mkdir -p ${EXECUTABLE_OUTPUT_PATH}/testQuadrature
#    )
