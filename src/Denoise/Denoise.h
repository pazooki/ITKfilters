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
//Filters
#include <itkConnectedThresholdImageFilter.h>
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include <itkBinaryBallStructuringElement.h>
#include <itkFlatStructuringElement.h>
#include <itkCurvatureAnisotropicDiffusionImageFilter.h>
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include <itkOtsuThresholdImageFilter.h>
#include <itkHuangThresholdImageFilter.h>
#include <itkYenThresholdImageFilter.h>
#include <itkShanbhagThresholdImageFilter.h>
class Denoise
{
public:
    // ITK Images typedefs;
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

    /**********************************************************/
    /*******Region Growth********/
    /**********************************************************/
    typedef itk::ConnectedThresholdImageFilter<RealImageType, RealImageType > ConnectedFilterType;
    typedef ConnectedFilterType::Pointer ConnectedFilterTypeP;

    /**********************************************************/
    /*******Anistropic Denoise********/
    /**********************************************************/
    typedef itk::GradientAnisotropicDiffusionImageFilter<
        RealImageType, RealImageType> AnisotropicFilterGradientType ;
    typedef AnisotropicFilterGradientType::Pointer AnisotropicFilterGradientTypeP;

    typedef itk::CurvatureAnisotropicDiffusionImageFilter<
        RealImageType, RealImageType> AnisotropicFilterCurvatureType ;
    typedef AnisotropicFilterCurvatureType::Pointer AnisotropicFilterCurvatureTypeP;

    /**********************************************************/
    /*******Morphological Operations********/
    /**********************************************************/
    // Ball
    typedef itk::BinaryBallStructuringElement< unsigned char, 2  > BallStructuringElementType;
    typedef itk::GrayscaleMorphologicalOpeningImageFilter<
        RealImageType , RealImageType, BallStructuringElementType >  BallOpeningFilterType;
    typedef BallOpeningFilterType::Pointer BallOpeningFilterTypeP;

    typedef itk::GrayscaleMorphologicalClosingImageFilter<
        RealImageType , RealImageType, BallStructuringElementType >  BallClosingFilterType;
    typedef BallClosingFilterType::Pointer BallClosingFilterTypeP;
    // Generic FlatStructure
    typedef itk::FlatStructuringElement< 2 > FlatStructuringElementType;
    typedef itk::GrayscaleMorphologicalOpeningImageFilter<
        RealImageType , RealImageType, FlatStructuringElementType >  FlatOpeningFilterType;
    typedef FlatOpeningFilterType::Pointer FlatOpeningFilterTypeP;
    typedef itk::GrayscaleMorphologicalClosingImageFilter<
        RealImageType , RealImageType, FlatStructuringElementType >  FlatClosingFilterType;
    typedef FlatClosingFilterType::Pointer FlatClosingFilterTypeP;

    /**********************************************************/
    /*******Binary Filters********/
    /**********************************************************/
    typedef itk::OtsuThresholdImageFilter<Denoise::RealImageType, Denoise::InputImageType> OtsuType;
    typedef OtsuType::Pointer OtsuTypeP;
    typedef itk::HuangThresholdImageFilter<Denoise::RealImageType, Denoise::InputImageType> HuangType;
    typedef HuangType::Pointer HuangTypeP;
    typedef itk::YenThresholdImageFilter<Denoise::RealImageType, Denoise::InputImageType> YenType;
    typedef YenType::Pointer YenTypeP;
    typedef itk::ShanbhagThresholdImageFilter<Denoise::RealImageType, Denoise::InputImageType> ShanbhagType;
    typedef ShanbhagType::Pointer ShanbhagTypeP;

public:
    Denoise() = default;
    // Denoise(const std::string &imgName);
    virtual ~Denoise(){};
    InputTypeP Read(const std::string &imgName);
    void Write(RealTypeP &input, std::string &imgName);
    void Write(InputTypeP &input, std::string &imgName);
    void Write(RealImageType* input, std::string &imgName);
    void Write(InputImageType* input, std::string &imgName);

    AnisotropicFilterGradientTypeP AnisotropicFilterGradient(RealTypeP img,
      const unsigned int numberOfIterations = 5, const double timeStep = 0.125 /* Recomm for 2D */ ,
      const double conductance = 2);
    AnisotropicFilterGradientTypeP AnisotropicFilterGradient(InputTypeP img,
      const unsigned int numberOfIterations = 5, const double timeStep = 0.125 /* Recomm for 2D */ ,
      const double conductance = 2 );
    AnisotropicFilterCurvatureTypeP AnisotropicFilterCurvature(RealTypeP img,
      const unsigned int numberOfIterations = 5, const double timeStep = 0.125 /* Recomm for 2D */ ,
      const double conductance = 2);
    AnisotropicFilterCurvatureTypeP AnisotropicFilterCurvature(InputTypeP img,
      const unsigned int numberOfIterations = 5, const double timeStep = 0.125 /* Recomm for 2D */ ,
      const double conductance = 2 );

    ConnectedFilterTypeP RegionGrowth(RealTypeP img,
            unsigned int lowThreshold, unsigned int highThreshold);
    ConnectedFilterTypeP RegionGrowth(InputTypeP img,
            unsigned int lowThreshold, unsigned int highThreshold);
    BallOpeningFilterTypeP MorphologicalOpening(RealTypeP img, int radius);
    BallOpeningFilterTypeP MorphologicalOpening(InputTypeP img, int radius);
    BallClosingFilterTypeP MorphologicalClosing(RealTypeP img, int radius);
    BallClosingFilterTypeP MorphologicalClosing(InputTypeP img, int radius);
    FlatOpeningFilterTypeP MorphologicalOpening(RealTypeP img, FlatStructuringElementType& structure );
    FlatOpeningFilterTypeP MorphologicalOpening(InputTypeP img, FlatStructuringElementType& structure);
    FlatClosingFilterTypeP MorphologicalClosing(RealTypeP img, FlatStructuringElementType& structure);
    FlatClosingFilterTypeP MorphologicalClosing(InputTypeP img, FlatStructuringElementType& structure);
    OtsuTypeP BinaryOtsu(RealTypeP img);
    OtsuTypeP BinaryOtsu(InputTypeP img);
    HuangTypeP BinaryHuang(RealTypeP img);
    HuangTypeP BinaryHuang(InputTypeP img);
    YenTypeP BinaryYen(RealTypeP img);
    YenTypeP BinaryYen(InputTypeP img);
    ShanbhagTypeP BinaryShanbhag(RealTypeP img);
    ShanbhagTypeP BinaryShanbhag(InputTypeP img);
public:
   InputTypeP inputImg_;
   // AnisotropicFilterTypeP anisotropicFilter_;
   ConnectedFilterTypeP filterGrow_;

};
#endif
