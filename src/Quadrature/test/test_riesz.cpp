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
#include <itkRescaleIntensityImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkSquareImageFilter.h>
#include <itkCastImageFilter.h>
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
    const string img_file{"./fixtures/collagen_98x98x20.tiff"};
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

    // Pad to valid FFT size. THIS SHIT put negative index values! Use FFTW to work with odd sizes. Or use hack of PasteImageFilter instead
    // typedef itk::FFTPadImageFilter<OutImageType> PadType;
    // auto padFilter = PadType::New();
    // padFilter->SetInput(castFilter->GetOutput());
    // padFilter->Update();
    // padFilter->UpdateLargestPossibleRegion();

    // Riesz filter
    typedef itk::RieszImageFilter<OutImageType> RieszFilter;
    auto riesz = RieszFilter::New();
    // riesz->SetSigmaGaussianDerivative(1.0);
    // riesz->SetInput(padFilter->GetOutput());
    riesz->SetInput(castFilter->GetOutput());
    riesz->Update();
    auto realComponent = riesz->GetOutputReal();
    auto rieszComponents = riesz->GetOutputRieszComponents();
    auto rieszNorm  = riesz->GetOutputRieszNorm();
    if(VFLAG) visualize::VisualizeITKImage(rieszNorm);

    auto eigenImage = riesz->ComputeEigenVectorsMaximizingRieszComponents(3, rieszComponents);

    // Compute Local Amplitude (sqrt(Even^2 + RieszNorm^2)

    // typedef itk::SquareImageFilter<OutImageType,OutImageType> SquareImageFilterType;
    // auto squareFilterEven = SquareImageFilterType::New();
    // squareFilterEven->SetInput(realComponent);
    // squareFilterEven->Update();
    // auto squareFilterNorm = SquareImageFilterType::New();
    // squareFilterNorm->SetInput(rieszNorm);
    // squareFilterNorm->Update();
    //
    // typedef itk::AddImageFilter<OutImageType,OutImageType, OutImageType> AddImageFilterType;
    // auto addFilter = AddImageFilterType::New();
    // addFilter->SetInput1(squareFilterEven->GetOutput());
    // addFilter->SetInput2(squareFilterNorm->GetOutput());
    // addFilter->Update();
    //
    // typedef SqrtImageFilter<OutImageType,OutImageType> SqrtImageFilterType;
    // typename SqrtImageFilterType::Pointer sqrtFilter = SqrtImageFilterType::New();
    // sqrtFilter->SetInput(addFilter->GetOutput());
    // sqrtFilter->Update();

    // if(VFLAG) visualize::VisualizeITKImages(sqrtFilter->GetOutput(),reader->GetOutput());
    // --- end of Local Amplitude ---


    // --- LocalPhase in unitary direction ---
    RieszFilter::DirectionType unitary_direction ;
    unitary_direction[0] = 1/sqrt(2.0);
    unitary_direction[1] = 0/sqrt(3.0);
    unitary_direction[2] = 0/sqrt(2.0);
    auto localPhaseDirection = riesz->ComputeLocalPhaseInDirection(unitary_direction, realComponent, rieszComponents);
    if(VFLAG) visualize::VisualizeITKImages(localPhaseDirection.GetPointer(), reader->GetOutput());
    // --- end of LocalPhase in unitary direction ---

    if(VFLAG) visualize::VisualizeITKImages(realComponent, reader->GetOutput());
    typedef VectorIndexSelectionCastImageFilter<typename RieszFilter::RieszComponentsImageType, OutImageType> CastIndexType;
    typename CastIndexType::Pointer castIndex= CastIndexType::New();
    castIndex->SetInput(rieszComponents);
    castIndex->SetIndex(0);
    if(VFLAG) visualize::VisualizeITKImage(castIndex->GetOutput());
    castIndex->SetIndex(1);
    if(VFLAG) visualize::VisualizeITKImage(castIndex->GetOutput());
    castIndex->SetIndex(2);
    castIndex->Update();
    if(VFLAG) visualize::VisualizeITKImage(castIndex->GetOutput());



}
