#include "gtest/gtest.h"
#include <memory>
#include <string>
#include <cmath>
#include "prog_options_test.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkForwardWaveletFilterBankFFT.h"
#include "itkHeldWavelet.h"
#include "itkForwardFFTImageFilter.h"
#include "itkInverseFFTImageFilter.h"
#include "visualize_functions.h"
// #include "itkWaveletFT.h"
using namespace std;
using namespace itk;

// TEST(quadrature, InCollagen){
int main(int argc, char** argv){
    auto option_map = program_options(argc, argv);
    bool VFLAG = option_map["visualize"].as<bool>();
    bool DEBUG = option_map["debug"].as<bool>();
    const string img_file{"/home/phc/repository_local/ITKfilters/src/fixtures/collagen_98x98x20.tiff"};
    const unsigned int dimension = 3;
    using PixelType = double;
    using ImageType = itk::Image<PixelType, dimension>;
    using ReaderType = itk::ImageFileReader<ImageType>;
    auto reader = ReaderType::New();
    reader->SetFileName(img_file);
    reader->Update();
    reader->UpdateLargestPossibleRegion();

    typedef itk::ForwardFFTImageFilter<ImageType> FFTFilterType;
    typename FFTFilterType::Pointer fftFilter = FFTFilterType::New();
    fftFilter->SetInput(reader->GetOutput());
    typedef typename FFTFilterType::OutputImageType ComplexImageType;

    typedef itk::HeldWavelet<> WaveletFunctionType;
    // One Step Filter Bank
    typedef itk::ForwardWaveletFilterBankFFT<ComplexImageType, ComplexImageType, WaveletFunctionType> ForwardFilterBank;
    typename ForwardFilterBank::Pointer forwardFilter = ForwardFilterBank::New();
    forwardFilter->SetHighPassSubBands(2);
    forwardFilter->SetInput(fftFilter->GetOutput());
    forwardFilter->SetDebug(DEBUG);
    forwardFilter->Update();
    auto lowPassImage = forwardFilter->GetOutput(0);
    auto highPassImage = forwardFilter->GetOutput(1);
    cout << lowPassImage->GetLargestPossibleRegion() << endl;
    cout << highPassImage->GetLargestPossibleRegion() << endl;

    // Inverse Transform
    typedef itk::InverseFFTImageFilter<ComplexImageType, ImageType> InverseFFTFilterType;
    typename InverseFFTFilterType::Pointer inverseFilter = InverseFFTFilterType::New();
    inverseFilter->SetInput(lowPassImage);
    inverseFilter->Update();
    if(VFLAG) visualize::VisualizeITKImage(inverseFilter->GetOutput());
    inverseFilter->SetInput(highPassImage);
    inverseFilter->Update();
    if(VFLAG) visualize::VisualizeITKImage(inverseFilter->GetOutput());

}

