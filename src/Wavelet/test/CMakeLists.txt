SET( ITK_LIBRARIES ${ITK_LIBRARIES} PARENT_SCOPE )
include_directories(SYSTEM ${GMOCK_INCLUDE_DIR} ${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Common)
include_directories("..")
include_directories("../WaveletFunctions")
# add_library(lpoTest prog_options_test.h)
# target_link_libraries(lpoTest ${Boost_PROGRAM_OPTIONS_LIBRARY})
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
find_package(Qt5 COMPONENTS Widgets)

# test_wavelet
add_executable(testWaveletFilterBankFFT.exe testWaveletFilterBankFFT.cpp)
set_source_files_properties (testWaveletFilterBankFFT.exe PROPERTIES COMPILE_FLAGS "-Wno-sign-compare")
target_link_libraries(testWaveletFilterBankFFT.exe gtest_main ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(testWaveletFilterBankFFT.exe ${ITK_LIBRARIES})
target_link_libraries(testWaveletFilterBankFFT.exe ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_link_libraries(testWaveletFilterBankFFT.exe ${VTK_LIBRARIES})
target_link_libraries(testWaveletFilterBankFFT.exe Qt5::Widgets)
message( ${ITK_LIBRARIES})
# add_custom_command(
#     TARGET test_quadrature.exe PRE_BUILD
#     COMMAND mkdir -p ${EXECUTABLE_OUTPUT_PATH}/testQuadrature
#    )
