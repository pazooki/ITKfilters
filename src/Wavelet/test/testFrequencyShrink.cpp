#include "gtest/gtest.h"
#include <memory>
#include <string>
#include <cmath>
#include "prog_options_test.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkForwardFFTImageFilter.h"
#include "itkInverseFFTImageFilter.h"
#include "visualize_functions.h"
#include <itkComplexToRealImageFilter.h>
#include <itkFrequencyShrinkImageFilter.h>
#include <itkShrinkImageFilter.h>
using namespace std;
using namespace itk;

int main(int argc, char** argv){
    auto option_map = program_options(argc, argv);
    bool VFLAG = option_map["visualize"].as<bool>();
    // bool DEBUG = option_map["debug"].as<bool>();
    // unsigned int input_n = option_map["input_n"].as<int>();
    // unsigned int input_l = option_map["input_l"].as<int>();
    const string img_file{"/home/phc/repository_local/ITKfilters/src/fixtures/collagen_64x64x16.tiff"};
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

    // ShinkFrequency
    typedef itk::FrequencyShrinkImageFilter<ComplexImageType, ComplexImageType> ShrinkType;
    typename ShrinkType::Pointer shrinkFilter = ShrinkType::New();
    shrinkFilter->SetInput(fftFilter->GetOutput());
    shrinkFilter->SetShrinkFactors(2);
    shrinkFilter->Update();
    // std::cout << "Post ShrinkFilter" << std::endl;
    // std::cout << "Buffered Region" << shrinkFilter->GetOutput()->GetBufferedRegion() << std::endl;
    // std::cout << "Region" << shrinkFilter->GetOutput()->GetLargestPossibleRegion() << std::endl;

    // // InverseFFT
    typedef itk::InverseFFTImageFilter<ComplexImageType, ImageType> InverseFFTFilterType;
    typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();
    inverseFFT->SetInput(shrinkFilter->GetOutput());
    inverseFFT->Update();
    std::cout << "FrequencyShrinker" <<std::endl;
    if(VFLAG) visualize::VisualizeITKImage(inverseFFT->GetOutput());
    std::cout << "Original" <<std::endl;
    if(VFLAG) visualize::VisualizeITKImage(reader->GetOutput());
    //Compare with regular srhink filter.
    typedef itk::ShrinkImageFilter<ImageType, ImageType> RegularShrinkType;
    typename RegularShrinkType::Pointer regularShrinkFilter = RegularShrinkType::New();
    regularShrinkFilter->SetInput(reader->GetOutput());
    regularShrinkFilter->SetShrinkFactors(2);
    regularShrinkFilter->Update();
    std::cout << "Regular shrinker" <<std::endl;
    if(VFLAG) visualize::VisualizeITKImage(regularShrinkFilter->GetOutput());

    // VISUALIZE Complex to real
    typedef itk::ComplexToRealImageFilter<ComplexImageType, ImageType> ComplexToRealFilter;
    typename ComplexToRealFilter::Pointer complexToRealFilter = ComplexToRealFilter::New();
    complexToRealFilter->SetInput(fftFilter->GetOutput() );
    complexToRealFilter->Update();
    if(VFLAG) visualize::VisualizeITKImage(complexToRealFilter->GetOutput());
    typename ComplexToRealFilter::Pointer complexToRealFilterShrink = ComplexToRealFilter::New();
    complexToRealFilterShrink->SetInput(shrinkFilter->GetOutput() );
    complexToRealFilterShrink->Update();
    if(VFLAG) visualize::VisualizeITKImage(complexToRealFilterShrink->GetOutput());

}

