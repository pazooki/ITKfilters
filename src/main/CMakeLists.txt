add_executable(denoise.exe main.cpp)
# include_directories(${IPP_LIB_FOLDER})
target_link_libraries(denoise.exe ldenoise)
