#include "gtest/gtest.h"
#include <memory>
#include <string>
#include <cmath>
#include "prog_options_test.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkForwardWaveletFilterBankFFT.h"
#include "itkInverseWaveletFilterBankFFT.h"
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
    typedef itk::ForwardWaveletFilterBankFFT<ComplexImageType, ComplexImageType, WaveletFunctionType> ForwardWaveletFilterBankType;
    typename ForwardWaveletFilterBankType::Pointer forwardWavelet = ForwardWaveletFilterBankType::New();
    unsigned int high_sub_bands = 2;
    forwardWavelet->SetHighPassSubBands(high_sub_bands);
    forwardWavelet->SetInput(fftFilter->GetOutput());
    forwardWavelet->SetDebug(DEBUG);
    forwardWavelet->Update();
    auto lowPassImage = forwardWavelet->GetOutput(0);
    auto highPassImage = forwardWavelet->GetOutput(high_sub_bands);
    cout << lowPassImage->GetLargestPossibleRegion() << endl;
    cout << highPassImage->GetLargestPossibleRegion() << endl;
    auto bandImage = forwardWavelet->GetOutputSubBand(1);
    cout << bandImage->GetLargestPossibleRegion() << endl;

    // Inverse FFT Transform
    typedef itk::InverseFFTImageFilter<ComplexImageType, ImageType> InverseFFTFilterType;
    typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();
    inverseFFT->SetInput(lowPassImage);
    inverseFFT->Update();
    if(VFLAG) visualize::VisualizeITKImage(inverseFFT->GetOutput());
    inverseFFT->SetInput(highPassImage);
    inverseFFT->Update();
    if(VFLAG) visualize::VisualizeITKImage(inverseFFT->GetOutput());
    inverseFFT->SetInput(bandImage);
    inverseFFT->Update();
    if(VFLAG) visualize::VisualizeITKImage(inverseFFT->GetOutput());

    //Inverse Wavelet Transform
    // typedef itk::InverseWaveletFilterBankFFT<ComplexImageType, ComplexImageType, WaveletFunctionType> InverseWaveletFilterBankType;
    // typename InverseWaveletFilterBankType::Pointer inverseWavelet = InverseWaveletFilterBankType::New();
    // inverseWavelet->SetHighPassSubBands(high_sub_bands);
    // inverseWavelet->SetInputs(forwardWavelet->GetOutputs());
    // inverseWavelet->SetDebug(DEBUG);
    // inverseWavelet->Update();

}

