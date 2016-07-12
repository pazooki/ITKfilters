#include "gtest/gtest.h"
#include <memory>
#include <string>
#include <cmath>
#include "prog_options_test.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkWaveletFrequencyForward.h"
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
    forwardWavelet->SetDebug(DEBUG);
    forwardWavelet->Print(std::cout);
    forwardWavelet->Update();

    for (unsigned int nout = 0 ; nout < forwardWavelet->GetNumberOfOutputs() ; ++nout)
    {
      unsigned int lv, b;
      std::tie(lv,b) = forwardWavelet->OutputIndexToLevelBand(nout);
      std::cout <<"InputIndex: " << nout << " --> lv:" << lv << " b:" << b << std::endl;
    }

    /* test OutputIndexToLevelBand */
    unsigned int ne = 0;
    {
      unsigned int lv, b;
      std::tie(lv,b) = forwardWavelet->OutputIndexToLevelBand(0);
      std::cout <<"inputindex: " << 0 << "lv:" << lv << " b:" << b << std::endl;
      if (lv != levels || b != 0)
        ++ne;
    }
    {
      unsigned int lv, b;
      std::tie(lv,b) = forwardWavelet->OutputIndexToLevelBand(high_sub_bands);
      std::cout <<"inputindex: " << high_sub_bands << "lv:" << lv << " b:" << b << std::endl;
      if (lv != 1 || b != high_sub_bands)
        ++ne;
    }
    if (ne != 0) std::cout << "ERROR in OutputIndexToLevelBand" << std::endl;

    // Inverse FFT Transform (Multilevel)
    typedef itk::InverseFFTImageFilter<ComplexImageType, ImageType> InverseFFTFilterType;
    typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();
    for (unsigned int level = 0 ; level < levels ; ++level)
    {
        for (unsigned int band = 0 ; band < high_sub_bands ; ++band)
        {
            if (level == 0 && band == 0) // Low pass
            {
                unsigned int n_output = 0;
                std::cout << "OutputIndex : " << n_output <<std::endl;
                std::cout << "Level: " << level + 1 << " / " << forwardWavelet->GetLevels() << std::endl;
                std::cout << "Band: " << band  << " / " << forwardWavelet->GetHighPassSubBands() << std::endl;
                std::cout << "Largest Region: " << forwardWavelet->GetOutput(n_output)->GetLargestPossibleRegion() << std::endl;
                std::cout << "Origin: " << forwardWavelet->GetOutput(n_output)->GetOrigin() << std::endl;
                std::cout << "Spacing: " << forwardWavelet->GetOutput(n_output)->GetSpacing() << std::endl;

                inverseFFT->SetInput(forwardWavelet->GetOutput(n_output) );
                inverseFFT->Update();
                if(VFLAG) visualize::VisualizeITKImage(inverseFFT->GetOutput());
            }
            unsigned int n_output = 1 + level * forwardWavelet->GetHighPassSubBands() + band;
            std::cout << "OutputIndex : " << n_output <<std::endl;
            std::cout << "Level: " << level + 1 << " / " << forwardWavelet->GetLevels() << std::endl;
            std::cout << "Band: " << band + 1 << " / " << forwardWavelet->GetHighPassSubBands() << std::endl;
            std::cout << "Largest Region: " << forwardWavelet->GetOutput(n_output)->GetLargestPossibleRegion() << std::endl;
            std::cout << "Origin: " << forwardWavelet->GetOutput(n_output)->GetOrigin() << std::endl;
            std::cout << "Spacing: " << forwardWavelet->GetOutput(n_output)->GetSpacing() << std::endl;

            inverseFFT->SetInput(forwardWavelet->GetOutput(n_output) );
            inverseFFT->Update();
            if(VFLAG) visualize::VisualizeITKImage(inverseFFT->GetOutput());

        }
    }

    //Get real part of complex image for visualization
    typedef itk::ComplexToRealImageFilter<ComplexImageType, ImageType> ComplexToRealFilter;
    typename ComplexToRealFilter::Pointer complexToRealFilter = ComplexToRealFilter::New();
    for (unsigned int level = 0 ; level < levels ; ++level)
    {
        for (unsigned int band = 0 ; band < high_sub_bands ; ++band)
        {
            if (level == 0 && band == 0) // Low pass
            {
                unsigned int n_output = 0;
                std::cout << "OutputIndex : " << n_output <<std::endl;
                std::cout << "Level: " << level + 1 << " / " << forwardWavelet->GetLevels() << std::endl;
                std::cout << "Band: " << band  << " / " << forwardWavelet->GetHighPassSubBands() << std::endl;
                std::cout << "Largest Region: " << forwardWavelet->GetOutput(n_output)->GetLargestPossibleRegion() << std::endl;
                std::cout << "Origin: " << forwardWavelet->GetOutput(n_output)->GetOrigin() << std::endl;
                std::cout << "Spacing: " << forwardWavelet->GetOutput(n_output)->GetSpacing() << std::endl;

                complexToRealFilter->SetInput(forwardWavelet->GetOutput(n_output) );
                complexToRealFilter->Update();
                // if(VFLAG) visualize::VisualizeITKImage(complexToRealFilter->GetOutput());
            }
            unsigned int n_output = 1 + level * forwardWavelet->GetHighPassSubBands() + band;
            std::cout << "OutputIndex : " << n_output <<std::endl;
            std::cout << "Level: " << level + 1 << " / " << forwardWavelet->GetLevels() << std::endl;
            std::cout << "Band: " << band + 1 << " / " << forwardWavelet->GetHighPassSubBands() << std::endl;
            std::cout << "Largest Region: " << forwardWavelet->GetOutput(n_output)->GetLargestPossibleRegion() << std::endl;
            std::cout << "Origin: " << forwardWavelet->GetOutput(n_output)->GetOrigin() << std::endl;
            std::cout << "Spacing: " << forwardWavelet->GetOutput(n_output)->GetSpacing() << std::endl;

            complexToRealFilter->SetInput(forwardWavelet->GetOutput(n_output) );
            complexToRealFilter->Update();
            // if(VFLAG) visualize::VisualizeITKImage(complexToRealFilter->GetOutput());

        }
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
