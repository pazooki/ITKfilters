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
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkNumberToString.h"

#include <string>

// boost::program_options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
namespace po = boost::program_options;

// Visualize for dev/debug purposes. Set in cmake file. Requires VTK
#include "itkViewImage.h"

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
  const unsigned int& inputLevels,
  const unsigned int& inputBands,
  const bool applySoftThreshold,
  const bool visualize,
  const double thresholdNumOfSigmas = 2.0)
{
  const unsigned int Dimension = VDimension;
  std::cout << "Dimension in run: " << Dimension << std::endl;

  typedef double                             PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;

  itk::NumberToString< unsigned int > n2s;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImage );
  reader->Update();

  typedef itk::ZeroDCImageFilter< ImageType > ZeroDCFilterType;
  typename ZeroDCFilterType::Pointer zeroDCFilter = ZeroDCFilterType::New();
  zeroDCFilter->SetInput( reader->GetOutput() );
  zeroDCFilter->Update();

  // Perform FFT on input image.
  typedef itk::ForwardFFTImageFilter< typename ZeroDCFilterType::OutputImageType >
    FFTForwardFilterType;
  typename FFTForwardFilterType::Pointer fftForwardFilter = FFTForwardFilterType::New();
  fftForwardFilter->SetInput( zeroDCFilter->GetOutput() );
  fftForwardFilter->Update();
  typedef typename FFTForwardFilterType::OutputImageType ComplexImageType;

  typedef itk::InverseFFTImageFilter< ComplexImageType, ImageType > InverseFFTFilterType;

  // Forward Wavelet
  typedef TWaveletFunction WaveletFunctionType;
  typedef itk::WaveletFrequencyFilterBankGenerator< ComplexImageType, WaveletFunctionType >
    WaveletFilterBankType;
  typedef itk::WaveletFrequencyForward< ComplexImageType, ComplexImageType, WaveletFilterBankType > ForwardWaveletType;
  // typedef itk::WaveletFrequencyForwardUndecimated< ComplexImageType, ComplexImageType, WaveletFilterBankType > ForwardWaveletType;
  typename ForwardWaveletType::Pointer forwardWavelet = ForwardWaveletType::New();
  unsigned int highSubBands = inputBands;
  unsigned int levels = inputLevels;
  forwardWavelet->SetHighPassSubBands( highSubBands );
  forwardWavelet->SetLevels( levels );
  forwardWavelet->SetInput( fftForwardFilter->GetOutput() );
  forwardWavelet->Update();
  typename ForwardWaveletType::OutputsType analysisWavelets =
    forwardWavelet->GetOutputs();

  // Apply Monogenic signal to wavelet results
  typedef itk::MonogenicSignalFrequencyImageFilter< ComplexImageType >
    MonogenicSignalFrequencyFilterType;
  typedef typename MonogenicSignalFrequencyFilterType::OutputImageType
    VectorMonoOutputType;
  typedef itk::VectorInverseFFTImageFilter< VectorMonoOutputType >
    VectorInverseFFTType;
  typedef itk::PhaseAnalysisSoftThresholdImageFilter< typename VectorInverseFFTType::OutputImageType >
    PhaseAnalysisFilter;

  typename ForwardWaveletType::OutputsType modifiedWavelets;
  unsigned int numberOfOutputs = forwardWavelet->GetNumberOfOutputs();
  for ( unsigned int i = 0; i < forwardWavelet->GetNumberOfOutputs(); ++i )
    {
    std::cout << "Output #: " << i << " / " << numberOfOutputs - 1 << std::endl;
    // if( i == numberOfOutputs - 1 ) // TODO Held does not modify approx image, but it does not generate better results.
    //   {
    //   modifiedWavelets.push_back( analysisWavelets[i] );
    //   continue;
    //   }
    typename MonogenicSignalFrequencyFilterType::Pointer monoFilter =
      MonogenicSignalFrequencyFilterType::New();
    typename VectorInverseFFTType::Pointer vecInverseFFT =
      VectorInverseFFTType::New();
    typename PhaseAnalysisFilter::Pointer phaseAnalyzer =
      PhaseAnalysisFilter::New();
    typename FFTForwardFilterType::Pointer fftForwardPhaseFilter =
      FFTForwardFilterType::New();

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
    bool visualizeCoefficients = false;
    if ( visualizeCoefficients )
    {
      for ( unsigned int i = 0; i < forwardWavelet->GetNumberOfOutputs(); ++i )
      {
        typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();
        inverseFFT->SetInput(analysisWavelets[i]);
        inverseFFT->Update();
        itk::Testing::ViewImage( inverseFFT->GetOutput(), "WaveletCoef: output #" + n2s(i) );
        inverseFFT->SetInput(modifiedWavelets[i]);
        inverseFFT->Update();
        itk::Testing::ViewImage( inverseFFT->GetOutput(), "WaveletCoef. PhaseAnalyzed #" + n2s(i) );
      }
    }
  } //end visualize

  typedef itk::WaveletFrequencyInverse< ComplexImageType, ComplexImageType, WaveletFilterBankType > InverseWaveletType;
  // typedef itk::WaveletFrequencyInverseUndecimated< ComplexImageType, ComplexImageType, WaveletFilterBankType > InverseWaveletType;
  typename InverseWaveletType::Pointer inverseWavelet = InverseWaveletType::New();
  inverseWavelet->SetHighPassSubBands( highSubBands );
  inverseWavelet->SetLevels( levels );
  inverseWavelet->SetInputs( modifiedWavelets );
  // The coefficients are now phases, do not apply reconstrucction factors.
  inverseWavelet->ApplyReconstructionFactorsOff();
  inverseWavelet->Update();

  typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();
  inverseFFT->SetInput( inverseWavelet->GetOutput() );
  inverseFFT->Update();

  if(visualize){
    itk::Testing::ViewImage( reader->GetOutput(), "Input Image" );
    itk::Testing::ViewImage( inverseFFT->GetOutput(), "Inverse Wavelet" );
  }

  // Cast To Float for save as tiff.
  typedef itk::Image< float, Dimension >                   ImageFloatType;
  typedef itk::CastImageFilter< ImageType, ImageFloatType> CastFloatType;
  typename CastFloatType::Pointer caster = CastFloatType::New();
  caster->SetInput(inverseFFT->GetOutput());
  caster->Update();

  // typedef itk::ImageFileWriter< typename InverseFFTFilterType::OutputImageType > WriterType;
  typedef itk::ImageFileWriter< ImageFloatType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  std::string appendString = "_L" + n2s(inputLevels) + "_B" + n2s(inputBands) + "_S" + n2s(thresholdNumOfSigmas);
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
    ( "output,o", po::value<std::string>()->required(), "Output path. More images will be generated based on this name. Include extension (.nrrd recommended)." )
    ( "levels,l", po::value<unsigned int>()->required(), "Number of Levels for the wavelet decomposition." )
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
      return false;
      }
    po::notify ( vm );
  } catch ( const std::exception& e ) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  const std::string inputImage = vm["input"].as<std::string>();
  const std::string outputImage = vm["output"].as<std::string>();
  const unsigned int inputLevels = vm["levels"].as<unsigned int>();
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

  typedef double                                         PixelType;

  if ( dimension == 2 )
    {
    const unsigned int ImageDimension = 2;
    typedef itk::Point< PixelType, ImageDimension >        PointType;

    if(verbose)
      std::cout << "dimension:" << ImageDimension << std::endl;

    typedef itk::HeldIsotropicWavelet< PixelType, ImageDimension, PointType >
      HeldIsotropicWaveletType;
    typedef itk::VowIsotropicWavelet< PixelType, ImageDimension, PointType >
      VowIsotropicWaveletType;
    typedef itk::SimoncelliIsotropicWavelet< PixelType, ImageDimension, PointType >
      SimoncelliIsotropicWaveletType;
    typedef itk::ShannonIsotropicWavelet< PixelType, ImageDimension, PointType >
      ShannonIsotropicWaveletType;
    if ( waveletFunction == "Held" )
      {
      return runRieszWaveletPhaseAnalysisTest< ImageDimension, HeldIsotropicWaveletType >( inputImage,
        outputImage,
        inputLevels,
        inputBands,
        applySoftThreshold,
        visualize,
        thresholdNumOfSigmas);
      }
    else if ( waveletFunction == "Vow" )
      {
      return runRieszWaveletPhaseAnalysisTest< ImageDimension, VowIsotropicWaveletType >( inputImage,
        outputImage,
        inputLevels,
        inputBands,
        applySoftThreshold,
        visualize,
        thresholdNumOfSigmas);
      }
    else if ( waveletFunction == "Simoncelli" )
      {
      return runRieszWaveletPhaseAnalysisTest< ImageDimension, SimoncelliIsotropicWaveletType >( inputImage,
        outputImage,
        inputLevels,
        inputBands,
        applySoftThreshold,
        visualize,
        thresholdNumOfSigmas);
      }
    else if ( waveletFunction == "Shannon" )
      {
      return runRieszWaveletPhaseAnalysisTest< ImageDimension, ShannonIsotropicWaveletType >( inputImage,
        outputImage,
        inputLevels,
        inputBands,
        applySoftThreshold,
        visualize,
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
    typedef itk::Point< PixelType, ImageDimension >        PointType;
    if(verbose)
      std::cout << "dimension:" << ImageDimension << std::endl;

    typedef itk::HeldIsotropicWavelet< PixelType, ImageDimension, PointType >
      HeldIsotropicWaveletType;
    typedef itk::VowIsotropicWavelet< PixelType, ImageDimension, PointType >
      VowIsotropicWaveletType;
    typedef itk::SimoncelliIsotropicWavelet< PixelType, ImageDimension, PointType >
      SimoncelliIsotropicWaveletType;
    typedef itk::ShannonIsotropicWavelet< PixelType, ImageDimension, PointType >
      ShannonIsotropicWaveletType;
    if ( waveletFunction == "Held" )
      {
      return runRieszWaveletPhaseAnalysisTest< ImageDimension, HeldIsotropicWaveletType >( inputImage,
        outputImage,
        inputLevels,
        inputBands,
        applySoftThreshold,
        visualize,
        thresholdNumOfSigmas);
      }
    else if ( waveletFunction == "Vow" )
      {
      return runRieszWaveletPhaseAnalysisTest< ImageDimension, VowIsotropicWaveletType >( inputImage,
        outputImage,
        inputLevels,
        inputBands,
        applySoftThreshold,
        visualize,
        thresholdNumOfSigmas);
      }
    else if ( waveletFunction == "Simoncelli" )
      {
      return runRieszWaveletPhaseAnalysisTest< ImageDimension, SimoncelliIsotropicWaveletType >( inputImage,
        outputImage,
        inputLevels,
        inputBands,
        applySoftThreshold,
        visualize,
        thresholdNumOfSigmas);
      }
    else if ( waveletFunction == "Shannon" )
      {
      return runRieszWaveletPhaseAnalysisTest< ImageDimension, ShannonIsotropicWaveletType >( inputImage,
        outputImage,
        inputLevels,
        inputBands,
        applySoftThreshold,
        visualize,
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
