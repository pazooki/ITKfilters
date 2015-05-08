#include "Denoise.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCastImageFilter.h>
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include <itkCurvatureAnisotropicDiffusionImageFilter.h>
using namespace std;

Denoise::InputTypeP Denoise::Read(const string &inputName){
    // Reader
    typedef itk::ImageFileReader< InputImageType > ReaderType;
    auto reader = ReaderType::New();

    // SCIFIO type
    // typedef itk::SCIFIOImageIO  ImageIOType;
    // TIFF type
    // typedef itk::TIFFImageIO ImageIOType;
    // auto scifioIO = ImageIOType::New();
    // reader->SetImageIO(scifioIO);
    reader->SetFileName( inputName );

    try {
        reader->Update();
    } catch( itk::ExceptionObject & excp ) {
        std::cerr << "Problem encountered while reading image file : " <<
            inputName << std::endl;
        throw;
    }

    // typedef itk::CastImageFilter< InputImageType, InputImageType > castFilterType;
    // castFilterType::Pointer castFilter = castFilterType::New();
    // castFilter->SetInput(reader->GetOutput());
    // castFilter->Update();

    inputImg_ = reader->GetOutput();

    // QuickView viewer;
    // viewer.AddImage(inputImg_);
    // viewer.Visualize();
    // inputImg_->Register();
    return inputImg_;
}

Denoise::RealTypeP Denoise::AnisotropicFilter(Denoise::RealTypeP img,
        const unsigned int numberOfIterations, const double timeStep,
        const double conductance, Denoise::AnisotropicFilterID id)
{
    switch(id){
        case AnisotropicFilterID::GRADIENT : {
            using AnisotropicFilterType = itk::GradientAnisotropicDiffusionImageFilter<
                Denoise::RealImageType, Denoise::RealImageType>;
            auto filter = AnisotropicFilterType::New();
            filter->SetInput( img );
            // numberOfIterations = 5;
            // const double timeStep = 0.125;
            //  const double conductance = 5;
            filter->SetNumberOfIterations( numberOfIterations );
            filter->SetTimeStep( timeStep );
            filter->SetConductanceParameter( conductance );
            filter->Update();
            anisotropicFilter_ = filter->GetOutput();
            return anisotropicFilter_;
            break;
        }
        case AnisotropicFilterID::CURVATURE : {
            typedef itk::CurvatureAnisotropicDiffusionImageFilter<
                Denoise::RealImageType, Denoise::RealImageType> AnisotropicFilterType;
            auto filter = AnisotropicFilterType::New();
            filter->SetInput( img );
            // numberOfIterations = 5;
            // const double timeStep = 0.125;
            //  const double conductance = 5;
            filter->SetNumberOfIterations( numberOfIterations );
            filter->SetTimeStep( timeStep );
            filter->SetConductanceParameter( conductance );
            filter->Update();
            anisotropicFilter_ = filter->GetOutput();
            return anisotropicFilter_;
            break;
        }
        default: throw std::runtime_error("Unknown id in AnisotropicFilter");
    }
}

Denoise::RealTypeP Denoise::AnisotropicFilter(Denoise::InputTypeP img,
        const unsigned int numberOfIterations, const double timeStep,
        const double conductance, Denoise::AnisotropicFilterID id)
{

    typedef itk::CastImageFilter< Denoise::InputImageType, Denoise::RealImageType > CastFilterType;
    auto filter = CastFilterType::New();
    // InputImageType i = *img->Clone();
    filter->SetInput(img);
    filter->Update();
    return AnisotropicFilter(filter->GetOutput(),
                numberOfIterations,timeStep, conductance, id);
}

Denoise::ConnectedFilterType::Pointer Denoise::RegionGrowth(Denoise::RealTypeP img,
                    unsigned int lth, unsigned int hth)
{
    filterGrow_ = ConnectedFilterType::New();
    filterGrow_->SetInput(img);
    filterGrow_->SetLower( lth );
    filterGrow_->SetUpper( hth );
    filterGrow_->SetReplaceValue( 255 );
    return filterGrow_;
}

Denoise::ConnectedFilterType::Pointer Denoise::RegionGrowth(Denoise::InputTypeP img,
                    unsigned int lth, unsigned int hth)
{
    typedef itk::CastImageFilter< Denoise::InputImageType, Denoise::RealImageType > CastFilterType;
    auto filter = CastFilterType::New();
    // InputImageType i = *img->Clone();
    filter->SetInput(img);
    filter->Update();
    return RegionGrowth(filter->GetOutput(),
                lth,hth);
}
