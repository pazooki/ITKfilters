#include "gtest/gtest.h"
#include <memory>
#include <string>
#include <cmath>
#include "prog_options_test.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkWaveletFrequencyFilterBankInverseGenerator.h"
#include "itkHeldWavelet.h"
#include "itkForwardFFTImageFilter.h"
#include "itkInverseFFTImageFilter.h"
#include <itkComplexToRealImageFilter.h>
#include "visualize_functions.h"
using namespace std;
using namespace itk;

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
    fftFilter->Update();
    typedef typename FFTFilterType::OutputImageType ComplexImageType;

    typedef itk::HeldWavelet<> WaveletFunctionType;
    // One Step Filter Bank
    typedef itk::WaveletFrequencyFilterBankInverseGenerator< ComplexImageType, WaveletFunctionType> InverseWaveletFilterBankType;
    typename InverseWaveletFilterBankType::Pointer inverseFilterBank = InverseWaveletFilterBankType::New();
    unsigned int high_sub_bands = 5;
    inverseFilterBank->SetHighPassSubBands(high_sub_bands);
    inverseFilterBank->SetSize(fftFilter->GetOutput()->GetLargestPossibleRegion().GetSize());
    inverseFilterBank->SetDebug(DEBUG);
    inverseFilterBank->Update();
    inverseFilterBank->Print(std::cout);
    auto lowPassImage = inverseFilterBank->GetOutput(0);
    cout << lowPassImage->GetLargestPossibleRegion() << endl;
    auto bandImage = inverseFilterBank->GetOutputSubBand(1);
    cout << bandImage->GetLargestPossibleRegion() << endl;
    auto highPassImage = inverseFilterBank->GetOutput(high_sub_bands);
    cout << highPassImage->GetLargestPossibleRegion() << endl;

    //Get real part of complex image for visualization
    typedef itk::ComplexToRealImageFilter<ComplexImageType, ImageType> ComplexToRealFilter;
    typename ComplexToRealFilter::Pointer complexToRealFilter = ComplexToRealFilter::New();
    for (unsigned int i = 0 ; i < high_sub_bands + 1 ; ++i)
    {
        if(VFLAG) std::cout << "Band: " << i << " / " << inverseFilterBank->GetHighPassSubBands() << std::endl;
        complexToRealFilter->SetInput(inverseFilterBank->GetOutput(i) );
        complexToRealFilter->Update();
        if(VFLAG) visualize::VisualizeITKImage(complexToRealFilter->GetOutput());
    }

    // Inverse FFT Transform
    typedef itk::InverseFFTImageFilter<ComplexImageType, ImageType> InverseFFTFilterType;
    typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();
    for (unsigned int i = 0 ; i < high_sub_bands + 1 ; ++i)
    {
        if(VFLAG) std::cout << "Band: " << i << " / " << inverseFilterBank->GetHighPassSubBands() << std::endl;
        inverseFFT->SetInput(inverseFilterBank->GetOutput(i) );
        inverseFFT->Update();
        if(VFLAG) visualize::VisualizeITKImage(inverseFFT->GetOutput());
    }
}
