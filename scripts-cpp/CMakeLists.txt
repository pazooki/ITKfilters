cmake_minimum_required(VERSION 3.5)
cmake_policy(SET CMP0012 NEW) # TO avoid ipcv lib error

project(ITKfilters_scripts)

# Following: https://discourse.itk.org/t/find-all-available-modules-from-groupio/351/5
# Since ITK4.13, ITKImageIO
set(use_itk_modules
    ITKFastMarching
    # ITKAnisotropicSmoothing # For fastMarching denoise
    ITKCommon
    ITKConnectedComponents
    ITKIOImageBase
    ITKImageFilterBase
    ITKLabelMap
    ITKRegionGrowing
    ITKThresholding
    ITKVtkGlue
    ITKImageIO
    ITKImageFrequency # ITK 5.0 (ViewImage)
    )
set(optional_itk_modules
    # AnisotropicDiffusionLBR
    IsotropicWavelets
    )

# DEV: ITK does not support OPTIONAL_COMPONENTS
# It wil fail if not found anyways.
find_package(ITK 5.0 REQUIRED
    COMPONENTS ${use_itk_modules}
    OPTIONAL_COMPONENTS ${optional_itk_modules}
    )
include(${ITK_USE_FILE})
message(status "** ITK Libraries:  ${ITK_LIBRARIES}")
# message(status "** ITK Includes:  ${ITK_INCLUDE_DIRS}")

# For visualize image only (kind of optional for the future)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(Boost COMPONENTS
        program_options
        filesystem
        REQUIRED )
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

include_directories("./include")

add_executable(fastMarchingSegmentation fastMarchingSegmentation.cxx)
target_link_libraries(fastMarchingSegmentation ${Boost_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(regionGrowingSegmentation regionGrowingSegmentation.cxx)
target_link_libraries(regionGrowingSegmentation ${Boost_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(regionGrowingConfidenceConnectedSegmentation regionGrowingConfidenceConnectedSegmentation.cxx)
target_link_libraries(regionGrowingConfidenceConnectedSegmentation ${Boost_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(imageInfo imageInfo.cxx)
target_link_libraries(imageInfo ${Boost_LIBRARIES} ${ITK_LIBRARIES})

add_executable(convertToNRRD convertToNRRD.cxx)
target_link_libraries(convertToNRRD ${Boost_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(runViewImage runViewImage.cxx)
target_link_libraries(runViewImage ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(compareImages itkTestingComparisonImageFilterTest.cxx)
target_include_directories(compareImages PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(compareImages ${ITK_LIBRARIES})# ${VTK_LIBRARIES})

# if(ITK_IsotropicWavelets)
    message(STATUS "Adding IsotropicWavelets scripts")
    add_executable(rieszWaveletPhaseAnalysis rieszWaveletPhaseAnalysis.cxx)
    target_link_libraries(rieszWaveletPhaseAnalysis ${Boost_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

    add_executable(fftPadPositiveIndex fftPadPositiveIndex.cxx)
    target_link_libraries(fftPadPositiveIndex ${Boost_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

    add_executable(rotateRieszWithStructureTensorDirections rotateRieszWithStructureTensorDirections.cxx)
    target_link_libraries(rotateRieszWithStructureTensorDirections ${Boost_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})
# endif()

# if(ITK_AnisotropicDiffusionLBR)
    # message(STATUS "Adding AnisotropicDiffusionLBR script")
    # add_executable(denoiseCoherenceAnisotropicLBR denoiseCoherenceAnisotropicLBR.cxx)
    # target_link_libraries(denoiseCoherenceAnisotropicLBR ${ITK_LIBRARIES} ${VTK_LIBRARIES})
# endif()
