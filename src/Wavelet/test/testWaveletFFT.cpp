#include "gtest/gtest.h"
#include <memory>
#include <string>
#include <cmath>
#include "prog_options_test.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkForwardWaveletFFT.h"
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
    typedef itk::ForwardWaveletFFT<ComplexImageType, ComplexImageType, WaveletFilterBankType> ForwardWaveletType;
    typename ForwardWaveletType::Pointer forwardWavelet = ForwardWaveletType::New();
    unsigned int high_sub_bands = input_n;
    unsigned int levels = input_l;
    forwardWavelet->SetHighPassSubBands(high_sub_bands);
    forwardWavelet->SetLevels(levels);
    forwardWavelet->SetInput(fftFilter->GetOutput());
    forwardWavelet->SetDebug(DEBUG);
    forwardWavelet->Print(std::cout);
    forwardWavelet->Update();

    //Get real part of complex image for visualization
    typedef itk::ComplexToRealImageFilter<ComplexImageType, ImageType> ComplexToRealFilter;
    typename ComplexToRealFilter::Pointer complexToRealFilter = ComplexToRealFilter::New();
    for (unsigned int i = 0 ; i < high_sub_bands + 1 ; ++i)
    {
        std::cout << "Band: " << i << " / " << forwardWavelet->GetHighPassSubBands() << std::endl;
        std::cout << "Largest Region: " << forwardWavelet->GetOutput(i)->GetLargestPossibleRegion() << std::endl;

        complexToRealFilter->SetInput(forwardWavelet->GetOutput(i) );
        complexToRealFilter->Update();
        if(VFLAG) visualize::VisualizeITKImage(complexToRealFilter->GetOutput());
    }

    // Inverse FFT Transform
    typedef itk::InverseFFTImageFilter<ComplexImageType, ImageType> InverseFFTFilterType;
    typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();
    for (unsigned int i = 0 ; i < high_sub_bands + 1 ; ++i)
    {
        if(VFLAG) std::cout << "Band: " << i << " / " << forwardWavelet->GetHighPassSubBands() << std::endl;
        inverseFFT->SetInput(forwardWavelet->GetOutput(i) );
        inverseFFT->Update();
        if(VFLAG) visualize::VisualizeITKImage(inverseFFT->GetOutput());
    }

    // // Inverse Wavelet Transform
    // typedef itk::InverseWaveletFilterBankFFT<ComplexImageType, ComplexImageType, WaveletFunctionType> InverseWaveletFilterBankType;
    // typename InverseWaveletFilterBankType::Pointer inverseWavelet = InverseWaveletFilterBankType::New();
    // inverseWavelet->SetHighPassSubBands(high_sub_bands);
    // inverseWavelet->SetInputs(forwardWavelet->GetOutputs());
    // inverseWavelet->SetDebug(DEBUG);
    // inverseWavelet->Update();
    // // Visualize inverse
    // inverseFFT->SetInput(inverseWavelet->GetOutput());
    // inverseFFT->Update();
    // // Compare these two in itk_test
    // if(VFLAG) visualize::VisualizeITKImage(inverseFFT->GetOutput());
    // // if(VFLAG) visualize::VisualizeITKImage(reader->GetOutput());
}
