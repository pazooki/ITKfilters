#include "gtest/gtest.h"
#include <memory>
#include <string>
#include <cmath>
#include "prog_options_test.h"
#include "visualize_functions.h"
#include "itkRieszImageFilter.h"
#include "itkFFTPadImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include <itkComplexToRealImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkSquareImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkIntensityWindowingImageFilter.h>

#include "itkWaveletOperator.h"
#include "itkWaveletFilterBank.h"
#include "itkWaveletTransform.h"
// bool VFLAG;
// int main(int argc, char** argv) {
//     ::testing::InitGoogleTest(&argc, argv);
//     auto option_map = program_options(argc, argv);
//     VFLAG = option_map["visualize"].as<bool>();
//     return RUN_ALL_TESTS();
// }
using namespace std;
using namespace itk;

// TEST(quadrature, InCollagen){
int main(int argc, char** argv){
    auto option_map = program_options(argc, argv);
    bool VFLAG = option_map["visualize"].as<bool>();
    const string img_file{"/home/phc/repository_local/ITKfilters/src/fixtures/collagen_98x98x20.tiff"};
    const unsigned int dimension = 3;
    using PixelType = unsigned int;
    using ImageType = itk::Image<PixelType, dimension>;
    using ReaderType = itk::ImageFileReader<ImageType>;
    auto reader = ReaderType::New();
    reader->SetFileName(img_file);
    reader->Update();
    reader->UpdateLargestPossibleRegion();

    using OutPixelType = double;
    using OutImageType = itk::Image<OutPixelType, dimension>;

    // Cast input image to real type.
    typedef itk::CastImageFilter< ImageType, OutImageType> CastFilterType;
    auto castFilter = CastFilterType::New();
    castFilter->SetInput(reader->GetOutput());
    castFilter->Update();

  /* Wavelet choice */
  const itk::Wavelet::Wavelet wavelet = itk::Wavelet::SYMLET8;

  /* Forward Transformation */
  typedef itk::WaveletOperator<wavelet, itk::Wavelet::FORWARD, PixelType, dimension>             WaveletOperator;
  typedef itk::WaveletFilterBank<ImageType, ImageType, WaveletOperator, itk::Wavelet::FORWARD>  ForwardFilterBank;
  typedef itk::WaveletTransform<ImageType, ImageType, ForwardFilterBank, itk::Wavelet::FORWARD> WaveletTransformType;

  unsigned int decimFactor = 1;
  typename WaveletTransformType::Pointer transform = WaveletTransformType::New();
  transform->SetInput(reader->GetOutput());
  // transform->SetNumberOfDecompositions(1);
  // transform->SetSubsampleImageFactor(decimFactor);
  transform->Update();

  // #<{(| Inverse Transformation |)}>#
  // typedef itk::WaveletOperator<wavelet, itk::Wavelet::INVERSE, PixelType, dimension>               InvWaveletOperator;
  // typedef itk::WaveletFilterBank<ImageType, ImageType, InvWaveletOperator, itk::Wavelet::INVERSE> InverseFilterBank;
  // typedef itk::WaveletTransform<ImageType, ImageType, InverseFilterBank, itk::Wavelet::INVERSE>   InvTransformType;
  //
  // typename InvTransformType::Pointer invTransform = InvTransformType::New();
  // invTransform->SetInput(transform->GetOutput());
  // invTransform->SetSubsampleImageFactor(decimFactor);
  // invTransform->Update();
    // Riesz filter
    // typedef itk::RieszImageFilter<OutImageType> RieszFilter;
    // auto riesz = RieszFilter::New();
    // // riesz->SetSigmaGaussianDerivative(1.0);
    // // riesz->SetInput(padFilter->GetOutput());
    // riesz->SetInput(castFilter->GetOutput());
    // riesz->Update();
    // // Get outputs
    // auto realComponent = riesz->GetOutputReal();
    // auto rieszComponents = riesz->GetOutputRieszComponents();
    // auto rieszNorm  = riesz->GetOutputRieszNorm();
    // auto fftForward = riesz->GetOutputFFT();
    //
    // typedef VectorIndexSelectionCastImageFilter<typename
    //     RieszFilter::RieszComponentsImageType, OutImageType> CastIndexType;
    // // Visualize Riesz Components
    // typename CastIndexType::Pointer castIndex= CastIndexType::New();
    // castIndex->SetInput(rieszComponents);
    // castIndex->SetIndex(0);
    // std::cout << "RieszComponent[0] (x)" <<  std::endl;
    // if(VFLAG) visualize::VisualizeITKImage(castIndex->GetOutput());
    // castIndex->SetIndex(1);
    // std::cout << "RieszComponent[1] (y)" <<  std::endl;
    // if(VFLAG) visualize::VisualizeITKImage(castIndex->GetOutput());
    // castIndex->SetIndex(2);
    // castIndex->Update();
    // std::cout << "RieszComponent[2] (z)" <<  std::endl;
    // if(VFLAG) visualize::VisualizeITKImage(castIndex->GetOutput());
}
