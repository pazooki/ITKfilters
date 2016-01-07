cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_ECLIPSE_VERSION=4.4)
# cmake ../src  -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_ECLIPSE_GENERATE_SOURCE_PROJECT=TRUE
cmake_policy(SET CMP0012 NEW) # TO avoid ipcv lib error
if (UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0  -Wall -std=c++11")
endif()

project(ITKfilters_p)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
find_package(ITK REQUIRED)
include_directories(${ITK_INCLUDE_DIRS})
include(${ITK_USE_FILE})
message(status "** ITK Libraries:  ${ITK_LIBRARIES}")
message(status "** ITK Includes:  ${ITK_INCLUDE_DIRS}")

find_package(Boost 1.54.0 COMPONENTS
        system filesystem
        program_options
        serialization
        REQUIRED )
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})


IF(ENABLE_PARALLEL)
    find_package(OpenMP REQUIRED)
    if(${OPENMP_FOUND})
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    else()
        message("OpenMP library not found")
    endif()

    add_definitions(-DENABLE_PARALLEL)
ENDIF(ENABLE_PARALLEL)
add_subdirectory(main)
add_subdirectory(Denoise)
add_subdirectory(Skeleton)
set(ENABLE_PARALLEL "ON" CACHE BOOL "Enable to use openMP")
################ GTEST ################
# Setup using ExternalProject + add_subdirectory.
#    from: http://crascit.com/2015/07/25/cmake-gtest/
# Download and unpack gtest at configure time
set(ENABLE_GOOGLE_TEST "ON" CACHE BOOL "Automatic download of gtest")
IF(ENABLE_GOOGLE_TEST)
    file(COPY "./fixtures" DESTINATION ${EXECUTABLE_OUTPUT_PATH} )
    configure_file(./fixtures/gtest-CMakeLists.txt gtest-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gtest-download )
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gtest-download )

    # Add gtest directly to our build
    add_subdirectory(${CMAKE_BINARY_DIR}/gtest-src
        ${CMAKE_BINARY_DIR}/gtest-build
        EXCLUDE_FROM_ALL )

    # Add the gtest include directory, since gtest
    # doesn't add that dependency to its gtest target
    target_include_directories(gtest INTERFACE
        "${CMAKE_BINARY_DIR}/gtest-src/include" )
ENDIF(ENABLE_GOOGLE_TEST)
# TEST
### Google TEST AND MOCK
# There must be an environment variable (in .bashrc or profile.d)
# export GMOCK_HOME=/path/to/gmock-1.7.0
# set(ENABLE_GOOGLE_TEST "ON" CACHE BOOL "Requires ENV variable GMOCK_HOME pointing to gmock parent folder")

#IF(ENABLE_GOOGLE_TEST)
#    file(COPY "./fixtures" DESTINATION ${EXECUTABLE_OUTPUT_PATH} )

#    set(GMOCK_HOME "$ENV{GMOCK_HOME}" CACHE PATH "Set this to gmock parent folder")
#    #FIND GMOCK
    # find_path(GTEST_INCLUDE_DIR
    #          gtest/gtest.h
    #          HINTS ${GMOCK_HOME}/gtest/include
    #          PATH_SUFFIXES gtest
    #          DOC "gtest include dir")
    # find_path(GMOCK_INCLUDE_DIR
    #          gmock/gmock.h
    #          HINTS ${GMOCK_HOME}/include
    #          PATH_SUFFIXES gmock
    #          DOC "gmock include dir")
    # include_directories(SYSTEM ${GMOCK_INCLUDE_DIR} ${GTEST_INCLUDE_DIR})

    # find_library(GTEST_LIBRARY
    #             NAMES gtest
    #                   libgtest
    #                   libgtest.a
    #             HINTS ${GMOCK_HOME}/mybuild/gtest
    #                   ${GMOCK_HOME}/gtest
    #                   ${GMOCK_HOME}/gtest/mybuild)
    # find_library(GTEST_LIBRARY_MAIN
    #             NAMES gtest_main
    #                   libgtest_main
    #                   libgtest.a
    #             HINTS ${GMOCK_HOME}/mybuild/gtest
    #                   ${GMOCK_HOME}/gtest
    #                   ${GMOCK_HOME}/gtest/mybuild)

    # find_library(GMOCK_LIBRARY
    #             NAMES gmock
    #                   libgmock
    #                   libgmock.a
    #             HINTS ${GMOCK_HOME}/mybuild
    #                   ${GMOCK_HOME})
    # find_library(GMOCK_LIBRARY_MAIN
    #             NAMES gmock_main
    #                   libgmock_main
    #                   libgmock.a
    #             HINTS ${GMOCK_HOME}/mybuild
    #                   ${GMOCK_HOME})
# ENDIF(ENABLE_GOOGLE_TEST)

