SET( ITK_LIBRARIES ${ITK_LIBRARIES} PARENT_SCOPE )
include_directories(SYSTEM ${GMOCK_INCLUDE_DIR} ${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Common)
include_directories("..")
include_directories("../IsotropicWaveletFrequencyFunctions")
# add_library(lpoTest prog_options_test.h)
# target_link_libraries(lpoTest ${Boost_PROGRAM_OPTIONS_LIBRARY})
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
find_package(Qt5 COMPONENTS Widgets)

# testsWaveletFilterBank
# Forward
add_executable(testWaveletFrequencyFilterBankForwardGenerator.exe testWaveletFrequencyFilterBankForwardGenerator.cpp)
set_source_files_properties (testWaveletFrequencyFilterBankForwardGenerator.exe PROPERTIES COMPILE_FLAGS "-Wno-sign-compare")
target_link_libraries(testWaveletFrequencyFilterBankForwardGenerator.exe gtest_main ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(testWaveletFrequencyFilterBankForwardGenerator.exe ${ITK_LIBRARIES})
target_link_libraries(testWaveletFrequencyFilterBankForwardGenerator.exe ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_link_libraries(testWaveletFrequencyFilterBankForwardGenerator.exe ${VTK_LIBRARIES})
target_link_libraries(testWaveletFrequencyFilterBankForwardGenerator.exe Qt5::Widgets)

# Inverse
add_executable(testWaveletFrequencyFilterBankInverseGenerator.exe testWaveletFrequencyFilterBankInverseGenerator.cpp)
set_source_files_properties (testWaveletFrequencyFilterBankInverseGenerator.exe PROPERTIES COMPILE_FLAGS "-Wno-sign-compare")
target_link_libraries(testWaveletFrequencyFilterBankInverseGenerator.exe gtest_main ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(testWaveletFrequencyFilterBankInverseGenerator.exe ${ITK_LIBRARIES})
target_link_libraries(testWaveletFrequencyFilterBankInverseGenerator.exe ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_link_libraries(testWaveletFrequencyFilterBankInverseGenerator.exe ${VTK_LIBRARIES})
target_link_libraries(testWaveletFrequencyFilterBankInverseGenerator.exe Qt5::Widgets)

##############

# testWavelet
# add_executable(testWaveletFFT.exe testWaveletFFT.cpp)
# set_source_files_properties (testWaveletFFT.exe PROPERTIES COMPILE_FLAGS "-Wno-sign-compare")
# target_link_libraries(testWaveletFFT.exe gtest_main ${CMAKE_THREAD_LIBS_INIT})
# target_link_libraries(testWaveletFFT.exe ${ITK_LIBRARIES})
# target_link_libraries(testWaveletFFT.exe ${Boost_PROGRAM_OPTIONS_LIBRARY})
# target_link_libraries(testWaveletFFT.exe ${VTK_LIBRARIES})
# target_link_libraries(testWaveletFFT.exe Qt5::Widgets)
# add_custom_command(
#     TARGET test_quadrature.exe PRE_BUILD
#     COMMAND mkdir -p ${EXECUTABLE_OUTPUT_PATH}/testQuadrature
#    )
