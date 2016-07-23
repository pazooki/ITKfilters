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
#include <itkFrequencyExpandImageFilter.h>
#include <itkExpandImageFilter.h>
using namespace std;
using namespace itk;

int main(int argc, char** argv){
    auto option_map = program_options(argc, argv);
    bool VFLAG = option_map["visualize"].as<bool>();
    bool DEBUG = option_map["debug"].as<bool>();
    unsigned int input_n = option_map["input_n"].as<int>();
    unsigned int input_l = option_map["input_l"].as<int>();
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
    typedef itk::FrequencyExpandImageFilter<ComplexImageType, ComplexImageType> ExpandType;
    typename ExpandType::Pointer expandFilter = ExpandType::New();
    expandFilter->SetInput(fftFilter->GetOutput());
    expandFilter->SetExpandFactors(2);
    expandFilter->Update();
    // std::cout << "Post ExpandFilter" << std::endl;
    // std::cout << "Buffered Region" << expandFilter->GetOutput()->GetBufferedRegion() << std::endl;
    // std::cout << "Region" << expandFilter->GetOutput()->GetLargestPossibleRegion() << std::endl;

    // // InverseFFT
    typedef itk::InverseFFTImageFilter<ComplexImageType, ImageType> InverseFFTFilterType;
    typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();
    inverseFFT->SetInput(expandFilter->GetOutput());
    inverseFFT->Update();
    std::cout << "FrequencyExpander" <<std::endl;
    if(VFLAG) visualize::VisualizeITKImage(inverseFFT->GetOutput());
    std::cout << "Original" <<std::endl;
    if(VFLAG) visualize::VisualizeITKImage(reader->GetOutput());

    //Compare with regular expand filter.
    typedef itk::ExpandImageFilter<ImageType, ImageType> RegularExpandType;
    typename RegularExpandType::Pointer regularExpandFilter = RegularExpandType::New();
    regularExpandFilter->SetInput(reader->GetOutput());
    regularExpandFilter->SetExpandFactors(2);
    regularExpandFilter->Update();
    std::cout << "Regular expander" <<std::endl;
    if(VFLAG) visualize::VisualizeITKImage(regularExpandFilter->GetOutput());

    // VISUALIZE Complex to real
    typedef itk::ComplexToRealImageFilter<ComplexImageType, ImageType> ComplexToRealFilter;
    typename ComplexToRealFilter::Pointer complexToRealFilter = ComplexToRealFilter::New();
    complexToRealFilter->SetInput(fftFilter->GetOutput() );
    complexToRealFilter->Update();
    if(VFLAG) visualize::VisualizeITKImage(complexToRealFilter->GetOutput());
    typename ComplexToRealFilter::Pointer complexToRealFilterExpand = ComplexToRealFilter::New();
    complexToRealFilterExpand->SetInput(expandFilter->GetOutput() );
    complexToRealFilterExpand->Update();
    if(VFLAG) visualize::VisualizeITKImage(complexToRealFilterExpand->GetOutput());

}

