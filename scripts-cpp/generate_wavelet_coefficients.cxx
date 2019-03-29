#include <iostream>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkWaveletFrequencyForward.h"
#include "itkWaveletFrequencyForwardUndecimated.h"
#include "itkWaveletFrequencyFilterBankGenerator.h"
#include "itkHeldIsotropicWavelet.h"
#include "itkVowIsotropicWavelet.h"
#include "itkSimoncelliIsotropicWavelet.h"
#include "itkShannonIsotropicWavelet.h"
#include "itkForwardFFTImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"

namespace fs = boost::filesystem;
namespace po = boost::program_options;

template< unsigned int VImageDimension, typename TWaveletFunction, bool TUseUndecimated >
int run( const std::string& inputImage,
         const std::string& outputFolder,
         const unsigned int& inputLevels,
         const unsigned int& inputBands,
         const std::string& waveletFunction,
         const bool writeSpatial,
         const bool writeComplex,
         const bool writeRealPart,
         const bool writeImaginaryPart)
{
  constexpr unsigned int ImageDimension = VImageDimension;
  using PixelType = double;
  using ImageType = itk::Image< PixelType, ImageDimension >;
  // Read Image
  using ReaderType = itk::ImageFileReader< ImageType >;
  auto reader = ReaderType::New();
  reader->SetFileName(inputImage);

  // Perform FFT on input image.
  using FFTFilterType = itk::ForwardFFTImageFilter< ImageType >;
  auto fftFilter = FFTFilterType::New();
  fftFilter->SetInput( reader->GetOutput() );

  using ComplexImageType = typename FFTFilterType::OutputImageType;

  // Set the WaveletFunctionType and the WaveletFilterBank
  using WaveletFunctionType = TWaveletFunction;
  using WaveletFilterBankType = itk::WaveletFrequencyFilterBankGenerator< ComplexImageType, WaveletFunctionType >;
  using ForwardWaveletType = typename std::conditional<TUseUndecimated,
        itk::WaveletFrequencyForwardUndecimated< ComplexImageType, ComplexImageType, WaveletFilterBankType >,
        itk::WaveletFrequencyForward< ComplexImageType, ComplexImageType, WaveletFilterBankType >
          >::type;

  unsigned int highSubBands = inputBands;
  unsigned int levels = inputLevels;
  auto forwardWavelet = ForwardWaveletType::New();
  forwardWavelet->SetHighPassSubBands( highSubBands );
  forwardWavelet->SetLevels(levels);
  forwardWavelet->SetInput(fftFilter->GetOutput());
  forwardWavelet->Update();

  // Write all the wavelet coefficients
  fs::path inputImage_path(inputImage);
  fs::path outputImageBase_path(outputFolder);
  outputImageBase_path = outputImageBase_path / inputImage_path.stem();
  int forwardWaveletNOutputs = forwardWavelet->GetNumberOfOutputs();
  for(int linearIndex = 0; linearIndex < forwardWaveletNOutputs; ++linearIndex) {
    auto level_band_pair = forwardWavelet->OutputIndexToLevelBand(linearIndex);
    auto outputImageBasePerLevel_path = outputImageBase_path;
    outputImageBasePerLevel_path += TUseUndecimated ? "_Undecimated" : "_Decimated";
    outputImageBasePerLevel_path += "_W" + waveletFunction +
      "_l" + std::to_string(level_band_pair.first) +
      "_b" + std::to_string(level_band_pair.second);

    if(writeComplex)
    { // Complex image
      using WriterType = itk::ImageFileWriter<ComplexImageType>;
      auto writer = WriterType::New();
      writer->SetInput(forwardWavelet->GetOutput(linearIndex));
      auto complexOutput_path = outputImageBasePerLevel_path;
      complexOutput_path += inputImage_path.extension();
      writer->SetFileName(complexOutput_path.c_str());
      writer->Update();
    }
    if(writeRealPart)
    { // Real Part
      using ComplexToRealFilter = itk::ComplexToRealImageFilter<ComplexImageType, ImageType>;
      auto complexToReal = ComplexToRealFilter::New();
      complexToReal->SetInput(forwardWavelet->GetOutput(linearIndex));
      using WriterRealType = itk::ImageFileWriter<ImageType>;
      auto writeReal = WriterRealType::New();
      writeReal->SetInput(complexToReal->GetOutput());
      auto complexToReal_path = outputImageBasePerLevel_path;
      complexToReal_path += "_real";
      complexToReal_path += inputImage_path.extension();
      writeReal->SetFileName(complexToReal_path.c_str());
      writeReal->Update();
    }
    if(writeImaginaryPart)
    { // Imaginary Part
      using ComplexToImaginaryFilter = itk::ComplexToImaginaryImageFilter<ComplexImageType, ImageType>;
      auto complexToImaginary = ComplexToImaginaryFilter::New();
      complexToImaginary->SetInput(forwardWavelet->GetOutput(linearIndex));
      using WriterImaginaryType = itk::ImageFileWriter<ImageType>;
      auto writeImaginary = WriterImaginaryType::New();
      writeImaginary->SetInput(complexToImaginary->GetOutput());
      auto complexToImaginary_path = outputImageBasePerLevel_path;
      complexToImaginary_path += "_imaginary";
      complexToImaginary_path += inputImage_path.extension();
      writeImaginary->SetFileName(complexToImaginary_path.c_str());
      writeImaginary->Update();
    }
    if(writeSpatial)
    { // Inverse FFT (spatial domain)
      using InverseFFTFilterType = itk::InverseFFTImageFilter< ComplexImageType, ImageType >;
      auto inverseFFT = InverseFFTFilterType::New();
      inverseFFT->SetInput(forwardWavelet->GetOutput(linearIndex));
      using WriterRealType = itk::ImageFileWriter<ImageType>;
      auto writeReal = WriterRealType::New();
      writeReal->SetInput(inverseFFT->GetOutput());
      auto waveletCoefficientSpatial_path = outputImageBasePerLevel_path;
      waveletCoefficientSpatial_path += "_spatial";
      waveletCoefficientSpatial_path += inputImage_path.extension();
      writeReal->SetFileName(waveletCoefficientSpatial_path.c_str());
      writeReal->Update();
    }
  }

  std::cout << "Completed." << std::endl;
  return EXIT_SUCCESS;
}

template< unsigned int VImageDimension, typename TWaveletFunction>
int
run_decimated_or_undecimated(
    const bool useUndecimated,
    const std::string& filename,
    const std::string& outputFolder,
    const unsigned int& levels,
    const unsigned int& bands,
    const std::string& waveletFunction,
    const bool writeSpatial,
    const bool writeComplex,
    const bool writeRealPart,
    const bool writeImaginaryPart)
{
  if(useUndecimated) {
    constexpr bool TUseUndecimated = true;
    return run<VImageDimension, TWaveletFunction, TUseUndecimated>(filename, outputFolder, levels, bands,
          waveletFunction, writeSpatial, writeComplex, writeRealPart, writeImaginaryPart);
  } else {
    constexpr bool TUseUndecimated = false;
    return run<VImageDimension, TWaveletFunction, TUseUndecimated>(filename, outputFolder, levels, bands,
          waveletFunction, writeSpatial, writeComplex, writeRealPart, writeImaginaryPart);
  }
}

int main(int argc, char* const argv[]){

  /*-------------- Parse command line -----------------------------*/
  po::options_description opt_desc ( "Allowed options are: " );
  opt_desc.add_options()( "help,h", "display this message." );
  opt_desc.add_options()( "input,i", po::value<std::string>()->required(), "Input image file." );
  opt_desc.add_options()( "outputFolder,o", po::value<std::string>()->required(), "Folder to export the results.");
  opt_desc.add_options()( "dimension,d", po::value<int>()->required(), "Dimensionality of input image. 2 or 3.");
  opt_desc.add_options()( "levels,l", po::value<int>()->default_value(1), "Number of multiscale levels for the wavelet analysis.");
  opt_desc.add_options()( "bands,b", po::value<int>()->default_value(1), "Number of high frequency sub-bands for each level.");
  opt_desc.add_options()( "waveletFunction,w", po::value<std::string>()->default_value("Simoncelli"), "Mother isotropic wavelet, options: Held, Simoncelli, Vow, Shannon.");
  opt_desc.add_options()( "writeSpatial,s",  po::bool_switch()->default_value(false), "Write spatial (after IFFT) output" );
  opt_desc.add_options()( "writeComplex,c",  po::bool_switch()->default_value(false), "Write complex output" );
  opt_desc.add_options()( "writeRealPart,r",  po::bool_switch()->default_value(false), "Write the real part of the complex output" );
  opt_desc.add_options()( "writeImaginaryPart",  po::bool_switch()->default_value(false), "Write the imaginary part of the complex output" );

  opt_desc.add_options()( "useUndecimated,u", po::bool_switch()->default_value(false), "Use undecimated wavelet pyramid (default is decimated).");
  opt_desc.add_options()( "verbose,v",  po::bool_switch()->default_value(false), "verbose output" );
  opt_desc.add_options()( "visualize,t", po::bool_switch()->default_value(false), "Visualize thin result. Requires VISUALIZE option at build");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, opt_desc), vm);
    if (vm.count ( "help" ) || argc<=1 )
    {
      std::cout << "Basic usage:\n" << opt_desc << "\n";
      return EXIT_SUCCESS;
    }
    po::notify ( vm );
  } catch ( const std::exception& e ) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  //Parse options
  std::string filename = vm["input"].as<std::string>();
  int dim = vm["dimension"].as<int>();
  int levels = vm["levels"].as<int>();
  int bands = vm["bands"].as<int>();
  std::string waveletFunction = vm["waveletFunction"].as<std::string>();
  bool verbose = vm["verbose"].as<bool>();
  bool writeSpatial = vm["writeSpatial"].as<bool>();
  bool writeComplex = vm["writeComplex"].as<bool>();
  bool writeRealPart = vm["writeRealPart"].as<bool>();
  bool writeImaginaryPart = vm["writeImaginaryPart"].as<bool>();
  bool useUndecimated = vm["useUndecimated"].as<bool>();
  bool visualize = vm["visualize"].as<bool>();

  std::string outputFolder = vm["outputFolder"].as<std::string>();
  const fs::path output_folder_path{outputFolder};
  if(!fs::exists(output_folder_path)) {
      std::cerr << "output folder doesn't exist : " << output_folder_path.string() << std::endl;
      throw po::validation_error(po::validation_error::invalid_option_value, "output_folder_path");
  }
  if ( vm.count("waveletFunction") &&
     (!( waveletFunction == "Held" || waveletFunction == "Vow" ||
         waveletFunction == "Simoncelli" || waveletFunction == "Shannon" ))
     ) {
    throw po::validation_error(po::validation_error::invalid_option_value, "waveletFunction");
  }
  if ( vm.count("dimension") && (!( dim == 2 || dim == 3 ))
     ) {
    throw po::validation_error(po::validation_error::invalid_option_value, "dimension");
  }
  /*----------- command line--- End Parse -----------------------------*/

  using FunctionValueType = double;
  if(dim == 3) {
    const unsigned int ImageDimension = 3;
    using HeldIsotropicWaveletType =
      itk::HeldIsotropicWavelet< FunctionValueType, ImageDimension>;
    using VowIsotropicWaveletType =
      itk::VowIsotropicWavelet< FunctionValueType, ImageDimension>;
    using SimoncelliIsotropicWaveletType =
      itk::SimoncelliIsotropicWavelet< FunctionValueType, ImageDimension>;
    using ShannonIsotropicWaveletType =
      itk::ShannonIsotropicWavelet< FunctionValueType, ImageDimension>;

    if(waveletFunction == "Simoncelli")
      return run_decimated_or_undecimated<ImageDimension, SimoncelliIsotropicWaveletType>(useUndecimated, filename, outputFolder, levels, bands,
          waveletFunction, writeSpatial, writeComplex, writeRealPart, writeImaginaryPart);
    else if(waveletFunction == "Held")
      return run_decimated_or_undecimated<ImageDimension, HeldIsotropicWaveletType>(useUndecimated, filename, outputFolder, levels, bands,
          waveletFunction, writeSpatial, writeComplex, writeRealPart, writeImaginaryPart);
    else if(waveletFunction == "Vow")
      return run_decimated_or_undecimated<ImageDimension, VowIsotropicWaveletType>(useUndecimated, filename, outputFolder, levels, bands,
          waveletFunction, writeSpatial, writeComplex, writeRealPart, writeImaginaryPart);
    else if(waveletFunction == "Shannon")
      return run_decimated_or_undecimated<ImageDimension, ShannonIsotropicWaveletType>(useUndecimated, filename, outputFolder, levels, bands,
          waveletFunction, writeSpatial, writeComplex, writeRealPart, writeImaginaryPart);
    else {
      std::cerr << "Unknown wavelet function: " << waveletFunction << std::endl;
      return EXIT_FAILURE;
    }
  } else {
    const unsigned int ImageDimension = 2;
    using HeldIsotropicWaveletType =
      itk::HeldIsotropicWavelet< FunctionValueType, ImageDimension>;
    using VowIsotropicWaveletType =
      itk::VowIsotropicWavelet< FunctionValueType, ImageDimension>;
    using SimoncelliIsotropicWaveletType =
      itk::SimoncelliIsotropicWavelet< FunctionValueType, ImageDimension>;
    using ShannonIsotropicWaveletType =
      itk::ShannonIsotropicWavelet< FunctionValueType, ImageDimension>;

    if(waveletFunction == "Simoncelli")
      return run_decimated_or_undecimated<ImageDimension, SimoncelliIsotropicWaveletType>(useUndecimated, filename, outputFolder, levels, bands,
          waveletFunction, writeSpatial, writeComplex, writeRealPart, writeImaginaryPart);
    else if(waveletFunction == "Held")
      return run_decimated_or_undecimated<ImageDimension, HeldIsotropicWaveletType>(useUndecimated, filename, outputFolder, levels, bands,
          waveletFunction, writeSpatial, writeComplex, writeRealPart, writeImaginaryPart);
    else if(waveletFunction == "Vow")
      return run_decimated_or_undecimated<ImageDimension, VowIsotropicWaveletType>(useUndecimated, filename, outputFolder, levels, bands,
          waveletFunction, writeSpatial, writeComplex, writeRealPart, writeImaginaryPart);
    else if(waveletFunction == "Shannon")
      return run_decimated_or_undecimated<ImageDimension, ShannonIsotropicWaveletType>(useUndecimated, filename, outputFolder, levels, bands,
          waveletFunction, writeSpatial, writeComplex, writeRealPart, writeImaginaryPart);
    else {
      std::cerr << "Unknown wavelet function: " << waveletFunction << std::endl;
      return EXIT_FAILURE;
    }
  }

}
