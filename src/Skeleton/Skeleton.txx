#ifndef Skeleton_txx
#define Skeleton_txx

#include "Skeleton.h"
using namespace std;
using namespace skeleton;
template< typename ImageType >
typename itk::ImageFileReader<ImageType>::Pointer skeleton::Read(const std::string & input)
{

    typedef itk::ImageFileReader< ImageType > ReaderType;
    auto reader = ReaderType::New();
    reader->SetFileName( input );
    try {
       reader->Update();
    } catch( itk::ExceptionObject & excp ) {
        std::cerr << "Problem encountered while reading image file : " <<
            input << std::endl;
        throw;
    }
    return reader;
};
itk::CurvatureAnisotropicDiffusionImageFilter<
        RealImageType, RealImageType> AnisotropicFilterCurvatureType ;
;
// int main( int argc, char* argv[] )
// {
//     std::cout << " A main" << std::endl;
//     string s{"hola"};
//     auto a = skeleton::Read<itk::Image<unsigned short, 3>>(s);
// }
//   if( argc < 3 )
//     {
//     std::cerr << "3D Version: "<< std::endl;
//     std::cerr << "Hessian + Multi-scale (Antagi): "<< std::endl;
//     std::cerr << " Hessian Parameters (SecondRankTensor are set for vessel extraction)";
//     std::cerr << " Scales for multi-scale analysis: 1 sigma is pixel size.";
//     std::cerr << "+ Centerline Extraction 3D (Homman): "<< std::endl;
//     std::cerr << "Usage: "<< std::endl;
//     std::cerr << argv[0];
//     std::cerr << " <InputFileName> <OutputFileName>";
//     std::cerr << " [SigmaMinimum] [SigmaMaximum]";
//     std::cerr << " [NumberOfSigmaSteps]";
//     std::cerr << std::endl;
//     return EXIT_FAILURE;
//     }
//
//   const char * inputFileName = argv[1];
//   const char * outputFileName = argv[2];
//   double sigmaMinimum = 1.0;
//   if( argc > 3 )
//     {
//     sigmaMinimum = atof( argv[3] );
//     }
//   double sigmaMaximum = 10.0;
//   if( argc > 4 )
//     {
//     sigmaMaximum = atof( argv[4] );
//     }
//   unsigned int numberOfSigmaSteps = 10;
//   if( argc > 5 )
//     {
//     numberOfSigmaSteps = atoi( argv[5] );
//     }
//
//   const unsigned int Dimension = 3;
//
//   typedef float                              PixelType;
//   typedef itk::Image< PixelType, Dimension > ImageType;
//
//   typedef itk::ImageFileReader< ImageType >  ReaderType;
//   ReaderType::Pointer reader = ReaderType::New();
//   reader->SetFileName( inputFileName );
//
//   typedef itk::SymmetricSecondRankTensor< double, Dimension > HessianPixelType;
//   typedef itk::Image< HessianPixelType, Dimension >           HessianImageType;
//   typedef itk::HessianToObjectnessMeasureImageFilter< HessianImageType, ImageType >
//     ObjectnessFilterType;
//   ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();
//   objectnessFilter->SetBrightObject( false );
//   objectnessFilter->SetScaleObjectnessMeasure( false );
//   objectnessFilter->SetAlpha( 0.5 );
//   objectnessFilter->SetBeta( 1.0 );
//   objectnessFilter->SetGamma( 5.0 );
//
//   typedef itk::MultiScaleHessianBasedMeasureImageFilter< ImageType, HessianImageType, ImageType >
//     MultiScaleEnhancementFilterType;
//   MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter =
//     MultiScaleEnhancementFilterType::New();
//   multiScaleEnhancementFilter->SetInput( reader->GetOutput() );
//   multiScaleEnhancementFilter->SetHessianToMeasureFilter( objectnessFilter );
//   multiScaleEnhancementFilter->SetSigmaStepMethodToLogarithmic();
//   multiScaleEnhancementFilter->SetSigmaMinimum( sigmaMinimum );
//   multiScaleEnhancementFilter->SetSigmaMaximum( sigmaMaximum );
//   multiScaleEnhancementFilter->SetNumberOfSigmaSteps( numberOfSigmaSteps );
//
//   typedef itk::Image< unsigned char, Dimension > OutputImageType;
//   typedef itk::RescaleIntensityImageFilter< ImageType, OutputImageType >
//     RescaleFilterType;
//   RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
//   rescaleFilter->SetInput( multiScaleEnhancementFilter->GetOutput() );
//
//   // Go to Binary from Hessian (MultiScale)
//   // Get Skeleton from Binary
//   typedef itk::BinaryThinningImageFilter <ImageType, ImageType> BinaryThinningImageFilterType;
//   BinaryThinningImageFilterType::Pointer binaryThinningImageFilter = BinaryThinningImageFilterType::New();
//   binaryThinningImageFilter->SetInput(image);
//   binaryThinningImageFilter->Update();
//
//   typedef itk::ImageFileWriter< OutputImageType > WriterType;
//   WriterType::Pointer writer = WriterType::New();
//   writer->SetFileName( outputFileName );
//   writer->SetInput( rescaleFilter->GetOutput() );
//   try
//     {
//     writer->Update();
//     }
//   catch( itk::ExceptionObject & error )
//     {
//     std::cerr << "Error: " << error << std::endl;
//     return EXIT_FAILURE;
//     }
//
//   return EXIT_SUCCESS;
// }
#endif
