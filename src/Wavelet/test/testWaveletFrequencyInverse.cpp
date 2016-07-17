#include "gtest/gtest.h"
#include <memory>
#include <string>
#include <cmath>
#include "prog_options_test.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkWaveletFrequencyForward.h"
#include "itkWaveletFrequencyInverse.h"
// #include "itkInverseWaveletFFT.h"
#include "itkWaveletFrequencyFilterBankGenerator.h"
#include "itkHeldWavelet.h"
#include "itkForwardFFTImageFilter.h"
#include "itkInverseFFTImageFilter.h"
#include "visualize_functions.h"
#include <itkComplexToRealImageFilter.h>
// #include "itkWaveletFT.h"
using namespace std;
using namespace itk;

// TEST(quadrature, InCollagen){
int main(int argc, char** argv){
    auto option_map = program_options(argc, argv);
    bool VFLAG = option_map["visualize"].as<bool>();
    bool DEBUG = option_map["debug"].as<bool>();
    unsigned int input_n = option_map["input_n"].as<int>();
    unsigned int input_l = option_map["input_l"].as<int>();
    const string img_file{"/home/phc/repository_local/ITKfilters/src/fixtures/collagen_98x98x20.tiff"};
    const unsigned int dimension = 3;
    using PixelType = double;
    using ImageType = itk::Image<PixelType, dimension>;
    using ReaderType = itk::ImageFileReader<ImageType>;
    auto reader = ReaderType::New();
    reader->SetFileName(img_file);
    reader->Update();
    reader->UpdateLargestPossibleRegion();

    // Perform FFT on input image.
    typedef itk::ForwardFFTImageFilter<ImageType> FFTFilterType;
    typename FFTFilterType::Pointer fftFilter = FFTFilterType::New();
    fftFilter->SetInput(reader->GetOutput());
    typedef typename FFTFilterType::OutputImageType ComplexImageType;

    // Set the WaveletFunctionType and the WaveletFilterBank
    typedef itk::HeldWavelet<> WaveletFunctionType;
    typedef itk::WaveletFrequencyFilterBankGenerator< ComplexImageType, WaveletFunctionType> WaveletFilterBankType;
    typedef itk::WaveletFrequencyForward<ComplexImageType, ComplexImageType, WaveletFilterBankType> ForwardWaveletType;
    typename ForwardWaveletType::Pointer forwardWavelet = ForwardWaveletType::New();
    unsigned int high_sub_bands = input_n;
    unsigned int levels = input_l;
    forwardWavelet->SetHighPassSubBands(high_sub_bands);
    forwardWavelet->SetLevels(levels);
    forwardWavelet->SetInput(fftFilter->GetOutput());
    forwardWavelet->Update();

    // Inverse Wavelet Transform
    typedef itk::WaveletFrequencyInverse<ComplexImageType, ComplexImageType, WaveletFilterBankType> InverseWaveletType;
    typename InverseWaveletType::Pointer inverseWavelet = InverseWaveletType::New();
    inverseWavelet->SetHighPassSubBands(high_sub_bands);
    inverseWavelet->SetLevels(levels);
    // inverseWavelet->SetInputLowPass(forwardWavelet->GetOutputLowPass())
    // inverseWavelet->SetInputsHighPass(forwardWavelet->GetOutputsHighPassBands())
    inverseWavelet->SetInputs(forwardWavelet->GetOutputs());
    inverseWavelet->SetDebug(DEBUG);
    inverseWavelet->Update();
    // Visualize inverse
    typedef itk::InverseFFTImageFilter<ComplexImageType, ImageType> InverseFFTFilterType;
    typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();
    inverseFFT->SetInput(inverseWavelet->GetOutput());
    inverseFFT->Update();
    // Compare these two in itk_test
    if(VFLAG) visualize::VisualizeITKImage(inverseFFT->GetOutput());
    if(VFLAG) visualize::VisualizeITKImage(reader->GetOutput());
}
