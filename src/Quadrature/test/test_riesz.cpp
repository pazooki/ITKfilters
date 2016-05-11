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
#include <itkIntensityWindowingImageFilter.h>
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

    // Pad to valid FFT size. THIS SHIT put negative index values!
    // Use FFTW to work with odd sizes. Or use hack of PasteImageFilter instead
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
    std::cout << "RieszNorm |Rx^2 + Ry^2 + Rz^2| of rieszComponents" << std::endl;
    if(VFLAG) visualize::VisualizeITKImage(rieszNorm);

    // auto eigenImage = riesz->ComputeEigenVectorsMaximizingRieszComponents(3, 20.0, rieszComponents);

    auto eigenPair = riesz->ComputeEigenAnalysisMaximizingRieszComponents(3, 20.0f, rieszComponents);
    auto maxLocalRiesz = riesz->ComputeRieszComponentsWithMaximumResponse(eigenPair.first, rieszComponents);
    // if(VFLAG) visualize::VisualizeITKImages(maxLocalRiesz.GetPointer(), reader->GetOutput());
    typedef VectorIndexSelectionCastImageFilter<typename RieszFilter::RieszComponentsImageType, OutImageType> CastIndexType;
    typename CastIndexType::Pointer castIndexMax= CastIndexType::New();
    castIndexMax->SetInput(maxLocalRiesz);
    castIndexMax->SetIndex(0);
    std::cout << "RieszComponents in local direction that locally maximizes response" <<  std::endl;
    std::cout << "RieszComponent[0] (first eigenVector)" <<  std::endl;
    if(VFLAG) visualize::VisualizeITKImage(castIndexMax->GetOutput());
    castIndexMax->SetIndex(1);
    std::cout << "RieszComponent[1] (second eigenVector)" <<  std::endl;
    if(VFLAG) visualize::VisualizeITKImage(castIndexMax->GetOutput());
    castIndexMax->SetIndex(2);
    castIndexMax->Update();
    std::cout << "RieszComponent[2] (third eigenVector)" <<  std::endl;
    if(VFLAG) visualize::VisualizeITKImage(castIndexMax->GetOutput());
    std::cout << "RieszNorm |Rx^2 + Ry^2 + Rz^2| of eigenAnalysis" << std::endl;
    RieszFilter::DirectionType weights;
    weights[0] = 1.0;
    weights[1] = 1.0;
    weights[1] = 1.0;
    std::cout << "    Uniformly weighted (1,1,1)" << std::endl;
    auto normEigen = riesz->ComputeRieszWeightedNorm(maxLocalRiesz, weights);
    if(VFLAG) visualize::VisualizeITKImage(normEigen.GetPointer());
    std::cout << "    Weighted by local eigenValues" << std::endl;
    auto normWeightedByEigenValues = riesz->ComputeRieszWeightedNormByEigenValues(rieszComponents, eigenPair.second );
    // if(VFLAG) visualize::VisualizeITKImage(normWeightedByEigenValues.GetPointer());
    typedef itk::IntensityWindowingImageFilter< OutImageType, OutImageType > WindowFilterType;
    typename WindowFilterType::Pointer windowFilter = WindowFilterType::New();
    windowFilter->SetInput(normWeightedByEigenValues);
    windowFilter->SetWindowMinimum(0);
    windowFilter->SetWindowMaximum(500000);
    windowFilter->SetOutputMinimum(0);
    windowFilter->SetOutputMaximum(500000);
    windowFilter->Update();
    std::cout << "Visualize Windowed" << std::endl;
    if(VFLAG) visualize::VisualizeITKImage(windowFilter->GetOutput());
    std::cout << "Amplitude of eigen (no windowed)" << std::endl;
    auto amplitudeEigen = riesz->ComputeLocalAmplitude(realComponent ,normWeightedByEigenValues );
    std::cout << "Amplitude of eigen (windowed)" << std::endl;
    auto amplitudeEigenWin = riesz->ComputeLocalAmplitude(realComponent ,windowFilter->GetOutput());
    if(VFLAG) visualize::VisualizeITKImage(amplitudeEigenWin.GetPointer());

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
    std::cout << "Local Phase of RieszComponents in direction:" << unitary_direction << std::endl;
    if(VFLAG) visualize::VisualizeITKImages(localPhaseDirection.GetPointer(), reader->GetOutput());
    // --- end of LocalPhase in unitary direction ---

    // if(VFLAG) visualize::VisualizeITKImages(realComponent, reader->GetOutput());
    typename CastIndexType::Pointer castIndex= CastIndexType::New();
    castIndex->SetInput(rieszComponents);
    castIndex->SetIndex(0);
    std::cout << "RieszComponent[0] (x)" <<  std::endl;
    if(VFLAG) visualize::VisualizeITKImage(castIndex->GetOutput());
    castIndex->SetIndex(1);
    std::cout << "RieszComponent[1] (y)" <<  std::endl;
    if(VFLAG) visualize::VisualizeITKImage(castIndex->GetOutput());
    castIndex->SetIndex(2);
    castIndex->Update();
    std::cout << "RieszComponent[2] (z)" <<  std::endl;
    if(VFLAG) visualize::VisualizeITKImage(castIndex->GetOutput());
}
