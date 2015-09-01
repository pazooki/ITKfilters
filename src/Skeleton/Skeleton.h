#ifndef Skeleton_h
#define Skeleton_h
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkHessianToObjectnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBinaryThinningImageFilter.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include <string>

#include <itkCurvatureAnisotropicDiffusionImageFilter.h>
namespace skeleton {
    const unsigned int Dimension = 3;

    template<typename ImageType>
    typename itk::ImageFileReader<ImageType>::Pointer Read(const std::string & input);

    template<typename ImageType>
    typename itk::ImageFileReader<ImageType>::Pointer Read(const std::string & input);

    template<typename ImageType>
    typename itk::CurvatureAnisotropicDiffusionImageFilter<
        ImageType, ImageType>
        AnisotropicFilterCurvature(
              ImageType img,
              const unsigned int numberOfIterations = 5,
              const double timeStep = 0.125 /* Recomm for 2D */ ,
              const double conductance = 2
              );
} // End Namespace
#endif
#include "Skeleton.txx"
