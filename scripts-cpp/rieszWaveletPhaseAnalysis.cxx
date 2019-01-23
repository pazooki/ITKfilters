/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkForwardFFTImageFilter.h"
#include "itkInverseFFTImageFilter.h"
#include "itkWaveletUtilities.h"
#include "itkWaveletFrequencyForward.h"
#include "itkWaveletFrequencyInverse.h"
#include "itkWaveletFrequencyForwardUndecimated.h"
#include "itkWaveletFrequencyInverseUndecimated.h"
#include "itkWaveletFrequencyFilterBankGenerator.h"
#include "itkHeldIsotropicWavelet.h"
#include "itkVowIsotropicWavelet.h"
#include "itkSimoncelliIsotropicWavelet.h"
#include "itkShannonIsotropicWavelet.h"

#include "itkMonogenicSignalFrequencyImageFilter.h"
#include "itkVectorInverseFFTImageFilter.h"
#include "itkPhaseAnalysisSoftThresholdImageFilter.h"
#include "itkZeroDCImageFilter.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkTIFFImageIO.h"
#include "itkChangeInformationImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkNumberToString.h"

#include <string>

// boost::program_options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
namespace po = boost::program_options;
// boost::filesystem
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

// Visualize for dev/debug purposes. Set in cmake file. Requires VTK
#include "itkViewImage.h"
#include <itkComplexToRealImageFilter.h>
#include <itkComplexToPhaseImageFilter.h>

std::string
AppendToFilenameRiesz(const std::string& filename, const std::string & appendix)
{
  std::size_t foundDot = filename.find_last_of('.');
  return filename.substr( 0, foundDot ) + appendix + filename.substr( foundDot );
}

// 1. Wavelet analysis (forward) on input image.
// 2. Create a Monogenic Signal (from Riesz function ) on each wavelet output..
// 3. Do a PhaseAnalysis on each Monogenic Signal.
// 4. Wavelet reconstruction (inverse) using as coefficients the output of the PhaseAnalysis.
// Without applying reconstruction factors: ApplyReconstructionFactorOff()
// 5. The result of the reconstruction will be an image that uses phase information at each level/band for improving local structure information, and can also work as an equalizator of brightness.
template< unsigned int VDimension, typename TWaveletFunction >
int
runRieszWaveletPhaseAnalysisTest( const std::string& inputImage,
  const std::string & outputImage,
  const std::string & inputLevels,
  const unsigned int inputBands,
  const bool applySoftThreshold,
  const bool visualize,
  const std::string & waveletFunction,
  const double thresholdNumOfSigmas = 2.0)
{
  const unsigned int Dimension = VDimension;

  using PixelType = double;
  using ImageType = itk::Image< PixelType, Dimension >;
  using ReaderType = itk::ImageFileReader< ImageType >;

  itk::NumberToString< unsigned int > n2s;
  auto reader = ReaderType::New();
  reader->SetFileName( inputImage );
  // Defense against dodgy tiff images (wrong metadata)
  // Read incorrectly with SCIFIO module.
  /* itk::TIFFImageIO::Pointer tiffIO = itk::TIFFImageIO::New(); */
  /* if(tiffIO->CanReadFile(inputImage.c_str())) */
  /*   reader->SetImageIO( tiffIO ); */
  reader->Update();

  auto sizeOriginal = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  unsigned int scaleFactor = 2;
  unsigned int maxNumberOfLevels = itk::utils::ComputeMaxNumberOfLevels(sizeOriginal, scaleFactor);
  std::cout << "MaxNumberOfLevels: " << maxNumberOfLevels  << ". Max level recommended: " << maxNumberOfLevels - 3 << std::endl;

  using ChangeInformationFilterType = itk::ChangeInformationImageFilter< ImageType >;
  auto changeInfoFilter = ChangeInformationFilterType::New();
  changeInfoFilter->SetInput( reader->GetOutput() );
  changeInfoFilter->ChangeDirectionOn();
  typename ImageType::DirectionType directionIdentity;
  directionIdentity.SetIdentity();
  changeInfoFilter->SetOutputDirection(directionIdentity);
  changeInfoFilter->Update();

  using ZeroDCFilterType = itk::ZeroDCImageFilter< ImageType >;
  auto zeroDCFilter = ZeroDCFilterType::New();
  zeroDCFilter->SetInput( changeInfoFilter->GetOutput() );
  zeroDCFilter->Update();

  // Perform FFT on input image.
  using FFTForwardFilterType =
    itk::ForwardFFTImageFilter< typename ZeroDCFilterType::OutputImageType >;
  auto fftForwardFilter = FFTForwardFilterType::New();
  fftForwardFilter->SetInput( zeroDCFilter->GetOutput() );
  fftForwardFilter->Update();
  using ComplexImageType = typename FFTForwardFilterType::OutputImageType;

  using InverseFFTFilterType = itk::InverseFFTImageFilter< ComplexImageType, ImageType >;

  // Forward Wavelet
  using WaveletFunctionType = TWaveletFunction;
  using WaveletFilterBankType = itk::WaveletFrequencyFilterBankGenerator< ComplexImageType, WaveletFunctionType >;
  using ForwardWaveletType = itk::WaveletFrequencyForward< ComplexImageType, ComplexImageType, WaveletFilterBankType >;
  // using ForwardWaveletType = itk::WaveletFrequencyForwardUndecimated< ComplexImageType, ComplexImageType, WaveletFilterBankType >;
  auto forwardWavelet = ForwardWaveletType::New();

  unsigned int levels = 0;
  if(inputLevels == "max"){
    levels = ForwardWaveletType::ComputeMaxNumberOfLevels(reader->GetOutput()->GetLargestPossibleRegion().GetSize(), forwardWavelet->GetScaleFactor() );
    std::cout << "maxLevels = " << levels << std::endl;
  } else {
    levels = std::stoul(inputLevels, nullptr, 0);
  }
  unsigned int highSubBands = inputBands;
  forwardWavelet->SetHighPassSubBands( highSubBands );
  forwardWavelet->SetLevels( levels );
  forwardWavelet->SetInput( fftForwardFilter->GetOutput() );
  forwardWavelet->Update();
  auto analysisWavelets = forwardWavelet->GetOutputs();

  // Apply Monogenic signal to wavelet results
  using MonogenicSignalFrequencyFilterType = itk::MonogenicSignalFrequencyImageFilter< ComplexImageType >;
  using VectorMonoOutputType = typename MonogenicSignalFrequencyFilterType::OutputImageType;
  using VectorInverseFFTType = itk::VectorInverseFFTImageFilter< VectorMonoOutputType >;
  using PhaseAnalysisFilter = itk::PhaseAnalysisSoftThresholdImageFilter< typename VectorInverseFFTType::OutputImageType >;

  typename ForwardWaveletType::OutputsType modifiedWavelets;
  unsigned int numberOfOutputs = forwardWavelet->GetNumberOfOutputs();
  for ( unsigned int i = 0; i < forwardWavelet->GetNumberOfOutputs(); ++i )
    {
    std::cout << "Output #: " << i << " / " << numberOfOutputs - 1 << std::endl;
    // if( i == numberOfOutputs - 1 )
    //   {
    //   modifiedWavelets.push_back( analysisWavelets[i] );
    //   continue;
    //   }
    auto monoFilter = MonogenicSignalFrequencyFilterType::New();
    auto vecInverseFFT = VectorInverseFFTType::New();
    auto phaseAnalyzer = PhaseAnalysisFilter::New();
    auto fftForwardPhaseFilter = FFTForwardFilterType::New();

    // Generate a monogenic signal (vector valued)
    monoFilter->SetInput( analysisWavelets[i] );
    monoFilter->Update();

    vecInverseFFT->SetInput( monoFilter->GetOutput() );
    vecInverseFFT->Update();

    phaseAnalyzer->SetInput( vecInverseFFT->GetOutput() );
    phaseAnalyzer->SetApplySoftThreshold( applySoftThreshold );
    if (applySoftThreshold)
      {
      phaseAnalyzer->SetNumOfSigmas(thresholdNumOfSigmas);
      }
    phaseAnalyzer->Update();

    fftForwardPhaseFilter->SetInput( phaseAnalyzer->GetOutputCosPhase() );
    fftForwardPhaseFilter->Update();

    modifiedWavelets.push_back( fftForwardPhaseFilter->GetOutput() );
    modifiedWavelets.back()->DisconnectPipeline();
    }

  if(visualize) {
    // Visualize and compare modified wavelets coefficients (and approx image)
    bool visualizeSpatialCoefficients = false;
    if ( visualizeSpatialCoefficients )
    {
      for ( unsigned int i = 0; i < forwardWavelet->GetNumberOfOutputs(); ++i )
      {
        auto inverseFFT = InverseFFTFilterType::New();
        inverseFFT->SetInput(analysisWavelets[i]);
        inverseFFT->Update();
        itk::ViewImage<ImageType>::View( inverseFFT->GetOutput(), "WaveletCoef (spatial domain) output #" + n2s(i) );
        inverseFFT->SetInput(modifiedWavelets[i]);
        inverseFFT->Update();
        itk::ViewImage<ImageType>::View( inverseFFT->GetOutput(), "WaveletCoef (spatial domain) PhaseAnalyzed #" + n2s(i) );
      }
    }
    bool visualizeFrequencyCoefficients = false;
    if ( visualizeFrequencyCoefficients )
    {
      for ( unsigned int i = 0; i < forwardWavelet->GetNumberOfOutputs(); ++i )
      {
        using ComplexToRealFilter = itk::ComplexToRealImageFilter< ComplexImageType, ImageType >;
        auto complexToRealFilter = ComplexToRealFilter::New();
        complexToRealFilter->SetInput(analysisWavelets[i]);
        complexToRealFilter->Update();
        itk::ViewImage<ImageType>::View( complexToRealFilter->GetOutput(), "WaveletCoef: frequency-domain RealPart #" + n2s(i) );

        using ComplexToPhaseFilter = itk::ComplexToPhaseImageFilter< ComplexImageType, ImageType >;
        auto complexToPhaseFilter = ComplexToPhaseFilter::New();
        complexToPhaseFilter->SetInput(analysisWavelets[i]);
        complexToPhaseFilter->Update();
        itk::ViewImage<ImageType>::View(complexToPhaseFilter->GetOutput(), "WaveletCoef: frequency-domain Phase #" + n2s(i));
      }
    }
  } //end visualize

  using InverseWaveletType = itk::WaveletFrequencyInverse< ComplexImageType, ComplexImageType, WaveletFilterBankType >;
  // using InverseWaveletType = itk::WaveletFrequencyInverseUndecimated< ComplexImageType, ComplexImageType, WaveletFilterBankType >;
  auto inverseWavelet = InverseWaveletType::New();
  inverseWavelet->SetHighPassSubBands( highSubBands );
  inverseWavelet->SetLevels( levels );
  inverseWavelet->SetInputs( modifiedWavelets );
  // The coefficients are now phases, do not apply reconstrucction factors.
  inverseWavelet->ApplyReconstructionFactorsOff();
  inverseWavelet->Update();

  auto inverseFFT = InverseFFTFilterType::New();
  inverseFFT->SetInput( inverseWavelet->GetOutput() );
  inverseFFT->Update();

  if(visualize){
    itk::ViewImage<ImageType>::View( reader->GetOutput(), "Input Image" );
    itk::ViewImage<ImageType>::View( inverseFFT->GetOutput(), "Inverse Wavelet" );
  }

  // Cast To Float for save as tiff.
  using ImageFloatType = itk::Image< float, Dimension >;
  using CastFloatType = itk::CastImageFilter< ImageType, ImageFloatType>;
  auto caster = CastFloatType::New();
  caster->SetInput(inverseFFT->GetOutput());
  caster->Update();

  // typedef itk::ImageFileWriter< typename InverseFFTFilterType::OutputImageType > WriterType;
  using WriterType = itk::ImageFileWriter< ImageFloatType >;
  auto writer = WriterType::New();
  std::string appendString = "_W"  + waveletFunction +
                             "_L"  + n2s(levels) +
                             "_B"  + n2s(inputBands);
  if(applySoftThreshold)
    appendString += "_ApplyS" + n2s(thresholdNumOfSigmas);
  std::string outputFile = AppendToFilenameRiesz(outputImage, appendString);
  writer->SetFileName( outputFile );
  writer->SetInput( caster->GetOutput() );

  writer->Update();
  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  /*-------------- Parse command line -----------------------------*/
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>()->required(), "Input image." )
    ( "outputFolder,o", po::value<std::string>()->required(), "Output folder path. A number of images will be generated depending on levels and bands." )
    ( "outputExtension,e", po::value<std::string>()->default_value("nrrd"), "Output extension." )
    ( "levels,l", po::value<std::string>()->required(), "Number of Levels for the wavelet decomposition. Allowed: positive digit or max" )
    ( "bands,b", po::value<unsigned int>()->required(), "Number of Bands for the wavelet decomposition." )
    ( "wavelet,w", po::value<std::string>()->default_value("Held"), "Type of Wavelet: Valid: Held, Simoncelli, Vow, Shannon." )
    ( "dimension,d", po::value<unsigned int>()->required(), "Dimension of the image: 2 or 3" )
    ( "apply", po::bool_switch()->default_value(false), "Apply Soft Threshold.")
    ( "threshold_sigmas,s", po::value<double>()->default_value(2.0), "Sigmas for soft threshold." )
    ( "visualize,t", po::bool_switch()->default_value(false), "Visualize using vtk based viewer.")
    ( "verbose,v",  po::bool_switch()->default_value(false), "verbose output." );

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
    if (vm.count ( "help" ) || argc<=1 )
      {
      std::cout << "Basic usage:\n" << general_opt << "\n";
      return EXIT_SUCCESS;
      }
    po::notify ( vm );
  } catch ( const std::exception& e ) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  const std::string inputImage = vm["input"].as<std::string>();
  const std::string outputFolder = vm["outputFolder"].as<std::string>();
  const std::string outputExtension = vm["outputExtension"].as<std::string>();
  const fs::path inputImageStem_path = fs::path(inputImage).stem();
  const fs::path outputFolder_path = fs::absolute(fs::path(outputFolder));
  const std::string outputImage = (outputFolder_path / fs::path(inputImageStem_path.string() + "." + outputExtension)).string();
  const std::string inputLevels = vm["levels"].as<std::string>();
  const unsigned int inputBands = vm["bands"].as<unsigned int>();
  const std::string waveletFunction = vm["wavelet"].as<std::string>();
  const unsigned int dimension = vm["dimension"].as<unsigned int>();
  const bool applySoftThreshold = vm["apply"].as<bool>();
  const double thresholdNumOfSigmas = vm["threshold_sigmas"].as<double>();
  const bool verbose = vm["verbose"].as<bool>();
  const bool visualize = vm["visualize"].as<bool>();

  if ( vm.count("wavelet") &&
     (!( waveletFunction == "Held" || waveletFunction == "Vow" ||
         waveletFunction == "Simoncelli" || waveletFunction == "Shannon" ))
     )
     throw po::validation_error(po::validation_error::invalid_option_value, "wavelet");
  // END PARSE

  using PixelType = double;

  if ( dimension == 2 )
    {
    const unsigned int ImageDimension = 2;
    using PointType = itk::Point< PixelType, ImageDimension >;

    using HeldIsotropicWaveletType =
      itk::HeldIsotropicWavelet< PixelType, ImageDimension, PointType >;
    using VowIsotropicWaveletType =
      itk::VowIsotropicWavelet< PixelType, ImageDimension, PointType >;
    using SimoncelliIsotropicWaveletType =
      itk::SimoncelliIsotropicWavelet< PixelType, ImageDimension, PointType >;
    using ShannonIsotropicWaveletType =
      itk::ShannonIsotropicWavelet< PixelType, ImageDimension, PointType >;
    if ( waveletFunction == "Held" )
      {
      return runRieszWaveletPhaseAnalysisTest<
        ImageDimension, HeldIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
            applySoftThreshold,
            visualize,
            waveletFunction,
            thresholdNumOfSigmas);
      }
    else if ( waveletFunction == "Vow" )
      {
      return runRieszWaveletPhaseAnalysisTest<
        ImageDimension, VowIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
            applySoftThreshold,
            visualize,
            waveletFunction,
            thresholdNumOfSigmas);
      }
    else if ( waveletFunction == "Simoncelli" )
      {
      return runRieszWaveletPhaseAnalysisTest<
        ImageDimension, SimoncelliIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
            applySoftThreshold,
            visualize,
            waveletFunction,
            thresholdNumOfSigmas);
      }
    else if ( waveletFunction == "Shannon" )
      {
      return runRieszWaveletPhaseAnalysisTest<
        ImageDimension, ShannonIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
            applySoftThreshold,
            visualize,
            waveletFunction,
            thresholdNumOfSigmas);
      }
    else
      {
      std::cerr << "Test failed!" << std::endl;
      std::cerr << waveletFunction << " wavelet type not supported." << std::endl;
      return EXIT_FAILURE;
      }
    }
  else if ( dimension == 3 )
    {
    const unsigned int ImageDimension = 3;
    using PointType = itk::Point< PixelType, ImageDimension >;

    using HeldIsotropicWaveletType =
      itk::HeldIsotropicWavelet< PixelType, ImageDimension, PointType >;
    using VowIsotropicWaveletType =
      itk::VowIsotropicWavelet< PixelType, ImageDimension, PointType >;
    using SimoncelliIsotropicWaveletType =
      itk::SimoncelliIsotropicWavelet< PixelType, ImageDimension, PointType >;
    using ShannonIsotropicWaveletType =
      itk::ShannonIsotropicWavelet< PixelType, ImageDimension, PointType >;
    if ( waveletFunction == "Held" )
      {
      return runRieszWaveletPhaseAnalysisTest<
        ImageDimension, HeldIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
            applySoftThreshold,
            visualize,
            waveletFunction,
            thresholdNumOfSigmas);
      }
    else if ( waveletFunction == "Vow" )
      {
      return runRieszWaveletPhaseAnalysisTest<
        ImageDimension, VowIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
            applySoftThreshold,
            visualize,
            waveletFunction,
            thresholdNumOfSigmas);
      }
    else if ( waveletFunction == "Simoncelli" )
      {
      return runRieszWaveletPhaseAnalysisTest<
        ImageDimension, SimoncelliIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
            applySoftThreshold,
            visualize,
            waveletFunction,
            thresholdNumOfSigmas);
      }
    else if ( waveletFunction == "Shannon" )
      {
      return runRieszWaveletPhaseAnalysisTest<
        ImageDimension, ShannonIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
            applySoftThreshold,
            visualize,
            waveletFunction,
            thresholdNumOfSigmas);
      }
    else
      {
      std::cerr << "Test failed!" << std::endl;
      std::cerr << waveletFunction << " wavelet type not supported." << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    std::cerr << "Test failed!" << std::endl;
    std::cerr << "Error: only 2 or 3 dimensions allowed, " << dimension << " selected." << std::endl;
    return EXIT_FAILURE;
    }
}
