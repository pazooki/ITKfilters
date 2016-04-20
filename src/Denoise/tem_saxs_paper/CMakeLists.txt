include_directories(SYSTEM ${GMOCK_INCLUDE_DIR} ${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Common)
include_directories("..")

add_executable(tem_saxs_denoise.exe tem_saxs_paper.cpp)
set_source_files_properties (tem_saxs_denoise.exe PROPERTIES COMPILE_FLAGS "-Wno-sign-compare")
target_link_libraries(tem_saxs_denoise.exe  ldenoise)
# target_link_libraries(tem_saxs_denoise.exe  ${GMOCK_LIBRARY_MAIN} ${GTEST_LIBRARY_MAIN} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(tem_saxs_denoise.exe gtest_main ${CMAKE_THREAD_LIBS_INIT})

add_custom_command(
    TARGET tem_saxs_denoise.exe PRE_BUILD
    COMMAND mkdir -p ${EXECUTABLE_OUTPUT_PATH}/resultsTemSaxsPaper
    COMMAND ln -sf /home/phc/Dropbox/Shared-Geelong-Palmerston/pectin/Pectin1_acid/pectin_1_ice_Montage_1045_12K_8_bit.tif ${EXECUTABLE_OUTPUT_PATH}/resultsTemSaxsPaper/
    COMMAND ln -sf /home/phc/Dropbox/Shared-Geelong-Palmerston/Carrageenan/Carrageenan_K/Montage_832.tif ${EXECUTABLE_OUTPUT_PATH}/resultsTemSaxsPaper/
    COMMAND ln -sf /home/phc/Dropbox/Shared-Geelong-Palmerston/Carrageenan/Carrageenan_Na/Montage_851.tif ${EXECUTABLE_OUTPUT_PATH}/resultsTemSaxsPaper/
   )
