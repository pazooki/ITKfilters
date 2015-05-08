#ifndef DENOISE_H_
#define DENOISE_H_

#include "itkImage.h"
#include <boost/filesystem.hpp>
#include <string>
#include <vector>
#include <array>
#include <utility>
#include <algorithm>
#include <cmath>
#include <chrono>
//RegionGrowth
#include <itkConnectedThresholdImageFilter.h>

class Denoise
{
public:
    // ITK typedefs;
    const static unsigned int  Dimension = 2;
    typedef unsigned int       InputPixelType;
    typedef double             RealPixelType;
    typedef short int          ShortPixelType;
    typedef itk::Image< InputPixelType, Dimension> InputImageType;
    typedef itk::Image< RealPixelType, Dimension>   RealImageType;
    typedef itk::Image< ShortPixelType, Dimension>   ShortImageType;
    typedef itk::Image< std::complex<RealPixelType>, Dimension> ComplexImageType;
    typedef InputImageType::Pointer InputTypeP;
    typedef RealImageType::Pointer RealTypeP;
    typedef ShortImageType::Pointer ShortTypeP;
    typedef ComplexImageType::Pointer ComplexTypeP;
    typedef itk::ConnectedThresholdImageFilter<RealImageType, RealImageType > ConnectedFilterType;
public:
    Denoise() = default;
    // Denoise(const std::string &imgName);
    virtual ~Denoise(){};
    InputTypeP Read(const std::string &imgName);
    enum AnisotropicFilterID {GRADIENT = 0, CURVATURE = 1};
    RealTypeP AnisotropicFilter(RealTypeP img,
      const unsigned int numberOfIterations = 5, const double timeStep = 0.125 /* Recomm for 2D */ ,
      const double conductance = 2, AnisotropicFilterID id = AnisotropicFilterID::GRADIENT);

    RealTypeP AnisotropicFilter(InputTypeP img,
      const unsigned int numberOfIterations = 5, const double timeStep = 0.125 /* Recomm for 2D */ ,
      const double conductance = 2, AnisotropicFilterID id = AnisotropicFilterID::GRADIENT);

    ConnectedFilterType::Pointer RegionGrowth(RealTypeP img,
            unsigned int lowThreshold, unsigned int highThreshold);
    ConnectedFilterType::Pointer RegionGrowth(InputTypeP img,
            unsigned int lowThreshold, unsigned int highThreshold);
public:
   InputTypeP inputImg_;
   RealTypeP anisotropicFilter_;
   ConnectedFilterType::Pointer filterGrow_;

};
#endif
