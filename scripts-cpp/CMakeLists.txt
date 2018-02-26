cmake_minimum_required(VERSION 3.5)
cmake_policy(SET CMP0012 NEW) # TO avoid ipcv lib error

project(ITKfilters_scripts)

# Following: https://discourse.itk.org/t/find-all-available-modules-from-groupio/351/5
# Since ITK4.13, ITKImageIO
set(use_itk_modules
    IsotropicWavelets
    ITKFastMarching
    ITKRegionGrowing
    ITKAnisotropicSmoothing # For fastMarching denoise
    ITKVtkGlue
    ITKImageIO
    )
find_package(ITK COMPONENTS ${use_itk_modules} REQUIRED)
include(${ITK_USE_FILE})
message(status "** ITK Libraries:  ${ITK_LIBRARIES}")
message(status "** ITK Includes:  ${ITK_INCLUDE_DIRS}")

# For visualize image only (kind of optional for the future)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(Boost COMPONENTS
        program_options
        filesystem
        REQUIRED )
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

include_directories("./include")

add_executable(rieszWaveletPhaseAnalysis rieszWaveletPhaseAnalysis.cxx)
target_link_libraries(rieszWaveletPhaseAnalysis ${Boost_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(fastMarchingSegmentation fastMarchingSegmentation.cxx)
target_link_libraries(fastMarchingSegmentation ${Boost_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(regionGrowingSegmentation regionGrowingSegmentation.cxx)
target_link_libraries(regionGrowingSegmentation ${Boost_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(regionGrowingConfidenceConnectedSegmentation regionGrowingConfidenceConnectedSegmentation.cxx)
target_link_libraries(regionGrowingConfidenceConnectedSegmentation ${Boost_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(runViewImage runViewImage.cxx)
target_link_libraries(runViewImage ${ITK_LIBRARIES} ${VTK_LIBRARIES})
