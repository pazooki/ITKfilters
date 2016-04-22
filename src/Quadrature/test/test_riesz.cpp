#include "gtest/gtest.h"
#include <memory>
#include <string>
#include "prog_options_test.h"
#include "visualize_functions.h"
#include "itkRieszImageFilter.h"
#include "itkFFTPadImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include <itkComplexToRealImageFilter.h>

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

    // Pad to valid FFT size.
    typedef itk::FFTPadImageFilter<OutImageType> PadType;
    auto padFilter = PadType::New();
    padFilter->SetInput(castFilter->GetOutput());
    padFilter->Update();
    padFilter->UpdateLargestPossibleRegion();

    // Riesz filter
    typedef itk::RieszImageFilter<OutImageType> RieszFilter;
    auto riesz = RieszFilter::New();
    riesz->SetInput(padFilter->GetOutput());

    typedef typename RieszFilter::OutputImageType ComplexImageType;

    typedef itk::ComplexToRealImageFilter<ComplexImageType , OutImageType> CastFromComplexFilterType;
    auto castComplexFilter = CastFromComplexFilterType::New();
    castComplexFilter->SetInput(riesz->GetOutput());
    castComplexFilter->Update();
    if(VFLAG) visualize::VisualizeITKImages(reader->GetOutput(), castComplexFilter->GetOutput());


}
