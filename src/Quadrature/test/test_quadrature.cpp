#include "gtest/gtest.h"
#include <memory>
#include <string>
#include "prog_options_test.h"
#include "visualize_functions.h"
#include "itkAdaptiveFiltering3DImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkFFTPadImageFilter.h"

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
    const string img_file{"./fixtures/collagen_101x99x20.tiff"};
    const unsigned int dimension = 3;
    using PixelType = unsigned int;
    using ImageType = itk::Image<PixelType, dimension>;
    using ReaderType = itk::ImageFileReader<ImageType>;
    auto reader = ReaderType::New();
    reader->SetFileName(img_file);
    reader->Update();
    reader->UpdateLargestPossibleRegion();

    using OutPixelType = float;
    using OutImageType = itk::Image<OutPixelType, dimension>;

    // Cast input image to real type.
    typedef itk::CastImageFilter< ImageType, OutImageType> CastFilterType;
    auto castFilter = CastFilterType::New();
    castFilter->SetInput(reader->GetOutput());
    castFilter->Update();
    castFilter->UpdateLargestPossibleRegion();

    // Pad to valid FFT size.
    typedef itk::FFTPadImageFilter<OutImageType> PadType;
    auto padFilter = PadType::New();
    padFilter->SetInput(castFilter->GetOutput());
    padFilter->Update();
    padFilter->UpdateLargestPossibleRegion();

    using AdaptiveType = itk::AdaptiveFiltering3DImageFilter<OutImageType,OutImageType>;
    auto filter = AdaptiveType::New();
    filter->SetInput(padFilter->GetOutput());
    filter->Update();
    filter->UpdateLargestPossibleRegion();
    if(VFLAG) visualize::VisualizeITKImages(reader->GetOutput(), filter->GetOutput());


}
