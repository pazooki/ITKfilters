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

#include "itkRieszFrequencyFilterBankGenerator.h"
#include "itkStructureTensor.h"
#include "itkRieszRotationMatrix.h"
#include "itkCompensatedSummation.h"

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
#include <itkComplexToImaginaryImageFilter.h>
#include <itkComplexToPhaseImageFilter.h>
#include <itkComplexToModulusImageFilter.h>

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
runRotateRieszWithStructureTensor( const std::string& inputImage,
  const std::string & outputImage,
  const std::string & inputLevels,
  const unsigned int inputBands,
  const unsigned int rieszOrder,
  const bool visualize,
  const std::string & waveletFunction
  )
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
  // using ForwardWaveletType = itk::WaveletFrequencyForward< ComplexImageType, ComplexImageType, WaveletFilterBankType >;
  using ForwardWaveletType = itk::WaveletFrequencyForwardUndecimated< ComplexImageType, ComplexImageType, WaveletFilterBankType >;
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

  // Create a Riesz filter bank (complex) to couple with each level of the wavelet
  using RieszFilterBankType = itk::RieszFrequencyFilterBankGenerator< ComplexImageType >;
  using MultiplyFilterType = itk::MultiplyImageFilter< ComplexImageType >;
  typename ForwardWaveletType::OutputsType modifiedWavelets;
  unsigned int numberOfOutputs = forwardWavelet->GetNumberOfOutputs();
  unsigned int highFrequencyOutputIndex = 0;
  std::cout << "Input Riesz Order: " << rieszOrder << std::endl;
  // Compute the rotation matrix in the highest frequency level.
  // This rotation, will then be applied to all levels with RieszRotationMatrix
  // TODO: Check the differences of rotation matrices you can get (check Dicente paper).
  // To keep directionality, probably you want the rieszOrder == 1 to get the matrix,
  // and then convert it to a SteerableMatrix of whatever rieszOrder the user decide (via RieszRotatioMatrix)
  // SetRotatioMatrix with R from the StructureTensor analysis of order 1. Hint: Use GetRotationalMatrixFromOutputMatrix
  // Multiply the RieszFrequencyFilterBankGenerator (of order N) with the wavelet pyramid. (M images, where M depends on the rieszOrder)
  // Rotate, aka, apply, the SteerableMatrix to the M images generated at each level
  //
  // For texture features (1 value per input image): Compute the sum of the energy of the (M*levels) images.
  // Energy should be the square of all the pixel values (TODO: any normalization?)
  /************ Get rotation matrix from StructureTensor (NOTE: at each location/pixel!!) *************/
  using StructureTensorType = itk::StructureTensor< ImageType >;
  auto tensorFirstOrder = StructureTensorType::New();
    {
    // Compute the rotation matrix from the following index:
    unsigned int i = highFrequencyOutputIndex;
    std::cout << "Computing 1st order StructureTensor. Output #: " << i << " / " << numberOfOutputs - 1 << std::endl;
    auto filterBank = RieszFilterBankType::New();
    filterBank->SetOutputParametersFromImage(analysisWavelets[i]);
    filterBank->SetOrder(1);
    filterBank->Update();
    unsigned int numberOfRieszOutputs = filterBank->GetNumberOfOutputs();
    std::cout << "RieszOutputs for the 1st order: " << numberOfRieszOutputs << std::endl;
    std::vector< typename ComplexImageType::Pointer > rieszOutputs = filterBank->GetOutputs();
    std::vector< typename ComplexImageType::Pointer > rieszWavelets;
    std::vector< typename ImageType::Pointer > rieszWaveletsSpatial;
    for ( unsigned int rieszComp = 0; rieszComp < numberOfRieszOutputs; ++rieszComp )
      {
      // Multiply wavelet with riesz.
      auto multiplyWaveletRiesz = MultiplyFilterType::New();
      multiplyWaveletRiesz->SetInput1(analysisWavelets[i]);
      multiplyWaveletRiesz->SetInput2(rieszOutputs[rieszComp]);
      multiplyWaveletRiesz->Update();
      rieszWavelets.push_back(multiplyWaveletRiesz->GetOutput());
      auto inverseFFT = InverseFFTFilterType::New();
      inverseFFT->SetInput(rieszWavelets[rieszComp]);
      inverseFFT->Update();
      rieszWaveletsSpatial.push_back(inverseFFT->GetOutput());
      if(visualize)
        {
        bool visualizeRieszWavelets = true;
        if ( visualizeRieszWavelets )
          {
          itk::NumberToString< unsigned int > n2s;
          itk::ViewImage<ImageType>::View( inverseFFT->GetOutput(),
            "RieszWaveletCoef: output #" + n2s(i) + " RieszComp: " + n2s(rieszComp) );
          }
        bool visualizeRieszWaveletsInFrequency = false;
        if ( visualizeRieszWaveletsInFrequency )
          {
          itk::NumberToString< unsigned int > n2s;
          using ComplexToRealFilterType = itk::ComplexToRealImageFilter< ComplexImageType, ImageType >;
          using ComplexToImaginaryFilterType = itk::ComplexToImaginaryImageFilter< ComplexImageType, ImageType >;
          auto complexToReal = ComplexToRealFilterType::New();
          auto complexToImaginary = ComplexToImaginaryFilterType::New();
          complexToReal->SetInput(rieszWavelets[rieszComp]);
          complexToReal->Update();
          itk::ViewImage<ImageType>::View( complexToReal->GetOutput(),
            "REAL:RieszWaveletCoef: output #" + n2s(i) + " RieszComp: " + n2s(rieszComp) );
          complexToImaginary->SetInput(rieszWavelets[rieszComp]);
          complexToImaginary->Update();
          itk::ViewImage<ImageType>::View( complexToImaginary->GetOutput(),
            "IMAGINARY:RieszWaveletCoef: output #" + n2s(i) + " RieszComp: " + n2s(rieszComp) );
          }
        } // end visualize
      }

    tensorFirstOrder->SetInputs( rieszWaveletsSpatial );
    // TODO sigma should be a parameter of the script!
    // tensorFirstOrder->SetGaussianWindowRadius(3);
    tensorFirstOrder->Update();
    }
  // tensorFirstOrder stores the eigenvector (aka rotation) matrix at each pixel.
  // We will apply this matrix to each level of the pyramid. Please note that this is
  // different to what the monogenic analysis was doing.

  using SteerMatrix = itk::RieszRotationMatrix<PixelType, Dimension >;
  using SpatialRotationMatrix = typename SteerMatrix::SpatialRotationMatrixType;
  /************ Now use the rotation matrix (from 1st order) to compute the steerable matrix *************/
  /************ That can be applied to a Riesz filterbank of any order. ************/
  using CompensatedSummationType = itk::CompensatedSummation<long double>;
  CompensatedSummationType totalEnergy;
  std::vector<long double> energyPerOutput;

  for (unsigned int i = 0; i < forwardWavelet->GetNumberOfOutputs(); ++i)
    {
    CompensatedSummationType totalEnergyAtOutput;
    unsigned long long totalSamplesAtOutput(0);
    std::cout << "Output #: " << i << " / " << numberOfOutputs - 1 << std::endl;
    //Optimization: riesz for the first level is already pre-computed, but recomputing.
    auto filterBank = RieszFilterBankType::New();
    filterBank->SetOutputParametersFromImage(analysisWavelets[i]);
    filterBank->SetOrder(rieszOrder);
    filterBank->Update();
    unsigned int numberOfRieszOutputs = filterBank->GetNumberOfOutputs();
    std::vector< typename ComplexImageType::Pointer > rieszOutputs = filterBank->GetOutputs();
    std::vector< typename ComplexImageType::Pointer > rieszWavelets;
    std::vector< typename ImageType::Pointer > rieszWaveletsSpatial;
    for ( unsigned int rieszComp = 0; rieszComp < numberOfRieszOutputs; ++rieszComp )
      {
      // Multiply wavelet with riesz.
      auto multiplyWaveletRiesz = MultiplyFilterType::New();
      multiplyWaveletRiesz->SetInput1(analysisWavelets[i]);
      multiplyWaveletRiesz->SetInput2(rieszOutputs[rieszComp]);
      multiplyWaveletRiesz->Update();
      rieszWavelets.push_back(multiplyWaveletRiesz->GetOutput());
      auto inverseFFT = InverseFFTFilterType::New();
      inverseFFT->SetInput(rieszWavelets[rieszComp]);
      inverseFFT->Update();
      rieszWaveletsSpatial.push_back(inverseFFT->GetOutput());
      } // end rieszCompoenents loop
    // Use RieszRotationMatrix, set rotation to the pre-computed with the first order
    // And Multiply with the rieszWavelets transform.
    // This will give you another set of rieszWavelets' that is rotated.
    // in this rotated version, the first (or the last?) index is the one with max energy.
    // You can:
    //  - visualize stuff
    //  - compute the Energy = Sum_over_i{|pixel_value_at_i|^2}

    // XXX: Note. This approach require Undecimated pyramid (to re-use the rotation matrix for all levels)
    // Iterate over the tensorFirstOrder output and get the rotation matrices.
    // compute the steer matrix based on the rotation matrix, and multiply it with the rieszWavelets coefficients.
    // (you can do that in the frequency domain)
    // the rotatedRieszWavelets are the outputs you want, you are not interested in reconstruction,
    // so extract the texture feature (1 value) from it.
    // Compute the energy with itkStatisticsImageFilter.h (SumOfSquares)

    std::cout << "Apply each rotation matrix (computed at each pixel in the high frequency level) to the rest of the pyramid.  At Output: " <<  i << std::endl;
    using StructureTensorImageIterator = itk::ImageScanlineConstIterator<typename StructureTensorType::OutputImageType>;
    auto tensorIt = StructureTensorImageIterator(tensorFirstOrder->GetOutput(), tensorFirstOrder->GetOutput()->GetLargestPossibleRegion());
    // using RieszBankImageIterator = itk::ImageScanlineConstIterator<ComplexImageType>;
    using RieszBankImageIterator = itk::ImageScanlineConstIterator<ImageType>;
    std::vector<RieszBankImageIterator> rieszBankIterators;
    for (unsigned int rieszComp = 0; rieszComp < numberOfRieszOutputs; ++rieszComp)
      {
      // rieszBankIterators.emplace_back(RieszBankImageIterator(rieszWavelets[rieszComp], rieszWavelets[rieszComp]->GetLargestPossibleRegion()));
      rieszBankIterators.emplace_back(RieszBankImageIterator(rieszWaveletsSpatial[rieszComp], rieszWaveletsSpatial[rieszComp]->GetLargestPossibleRegion()));
      }
    while ( !tensorIt.IsAtEnd() )
      {
      while ( !tensorIt.IsAtEndOfLine() )
        {
        auto rotationMatrixVariableSize = tensorFirstOrder->GetRotationMatrixFromOutputMatrix(tensorIt.Get());
        SpatialRotationMatrix rotationMatrix;
        for (unsigned int row = 0; row < Dimension; ++row)
          {
          rotationMatrix.GetVnlMatrix().set_row(row, rotationMatrixVariableSize.GetVnlMatrix().get_row(row));
          }
        SteerMatrix steerMatrix;
        steerMatrix.SetOrder(rieszOrder);
        steerMatrix.SetSpatialRotationMatrix(rotationMatrix);
        steerMatrix.ComputeSteerableMatrix(); // Update()

        // std::vector<std::complex<PixelType>> inputRieszVector;
        std::vector<PixelType> inputRieszVector;
        for (unsigned int rieszComp = 0; rieszComp < numberOfRieszOutputs; ++rieszComp)
          {
          inputRieszVector.push_back(rieszBankIterators[rieszComp].Get() );
          }
        // auto rotatedRieszWavelets = steerMatrix.template MultiplyWithVector<std::complex<PixelType>>(inputRieszVector);
        auto rotatedRieszWavelets = steerMatrix.template MultiplyWithVector<PixelType>(inputRieszVector);
        for(auto & c : rotatedRieszWavelets)
          {
          // totalEnergyAtOutput += std::norm(c);
          totalEnergyAtOutput += c*c;
          ++totalSamplesAtOutput;
          }
        ++tensorIt;
        for (unsigned int rieszComp = 0; rieszComp < numberOfRieszOutputs; ++rieszComp)
          {
          ++rieszBankIterators[rieszComp];
          }
        }
      tensorIt.NextLine();
      for (unsigned int rieszComp = 0; rieszComp < numberOfRieszOutputs; ++rieszComp)
        {
        rieszBankIterators[rieszComp].NextLine();
        }
      }
    std::cout << "BEFORE: " << std::endl;
    energyPerOutput.push_back(totalEnergyAtOutput.GetSum() / totalSamplesAtOutput);
    std::cout << "energyPerOutput: " << energyPerOutput.back() << std::endl;
    totalEnergy+= energyPerOutput.back();
    // std::cout << "accumulatedEnergy: " << totalEnergy.GetSum() << std::endl;
    } // end forwardWavelet output


  long double totalEnergyD = totalEnergy.GetSum();
  std::cout << "ENERGY:" << std::endl;
  std::cout << "totalEnergy = " << totalEnergyD << std::endl;
  std::cout << "numberOfOutputs = " << numberOfOutputs << std::endl;
  std::cout << "energyPerOutput: " << std::endl;
  for (unsigned int energyIndex = 0; energyIndex < energyPerOutput.size(); ++energyIndex)
    {
    auto lv_band = forwardWavelet->OutputIndexToLevelBand(energyIndex);
    std::cout << "Output: " << energyIndex<< " L: " << lv_band.first << " B: " << lv_band.second << " EnergyPerOutput: " << energyPerOutput[energyIndex] << std::endl;
    }
  std::cout << "totalEnergy:Average = " << totalEnergyD / numberOfOutputs << std::endl;
  // typedef itk::ImageFileWriter< typename InverseFFTFilterType::OutputImageType > WriterType;
  // using WriterType = itk::ImageFileWriter< ImageFloatType >;
  // auto writer = WriterType::New();
  // std::string appendString = "_W"  + waveletFunction +
  //                            "_L"  + n2s(levels) +
  //                            "_B"  + n2s(inputBands);
  // std::string outputFile = AppendToFilenameRiesz(outputImage, appendString);
  // writer->SetFileName( outputFile );
  // writer->SetInput( caster->GetOutput() );
  //
  // writer->Update();
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
    ( "rieszOrder,r", po::value<unsigned int>()->default_value(1), "Order of the Riesz transform. 1 is equivalent to Monogenic signal." )
    ( "wavelet,w", po::value<std::string>()->default_value("Held"), "Type of Wavelet: Valid: Held, Simoncelli, Vow, Shannon." )
    ( "dimension,d", po::value<unsigned int>()->required(), "Dimension of the image: 2 or 3" )
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
  const unsigned int rieszOrder = vm["rieszOrder"].as<unsigned int>();
  const std::string waveletFunction = vm["wavelet"].as<std::string>();
  const unsigned int dimension = vm["dimension"].as<unsigned int>();
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
      return runRotateRieszWithStructureTensor<
        ImageDimension, HeldIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
			rieszOrder,
            visualize,
            waveletFunction
            );
      }
    else if ( waveletFunction == "Vow" )
      {
      return runRotateRieszWithStructureTensor<
        ImageDimension, VowIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
			rieszOrder,
            visualize,
            waveletFunction
            );
      }
    else if ( waveletFunction == "Simoncelli" )
      {
      return runRotateRieszWithStructureTensor<
        ImageDimension, SimoncelliIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
			rieszOrder,
            visualize,
            waveletFunction
            );
      }
    else if ( waveletFunction == "Shannon" )
      {
      return runRotateRieszWithStructureTensor<
        ImageDimension, ShannonIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
			rieszOrder,
            visualize,
            waveletFunction
            );
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
      return runRotateRieszWithStructureTensor<
        ImageDimension, HeldIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
			rieszOrder,
            visualize,
            waveletFunction
            );
      }
    else if ( waveletFunction == "Vow" )
      {
      return runRotateRieszWithStructureTensor<
        ImageDimension, VowIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
			rieszOrder,
            visualize,
            waveletFunction
            );
      }
    else if ( waveletFunction == "Simoncelli" )
      {
      return runRotateRieszWithStructureTensor<
        ImageDimension, SimoncelliIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
			rieszOrder,
            visualize,
            waveletFunction
            );
      }
    else if ( waveletFunction == "Shannon" )
      {
      return runRotateRieszWithStructureTensor<
        ImageDimension, ShannonIsotropicWaveletType >(
            inputImage,
            outputImage,
            inputLevels,
            inputBands,
			rieszOrder,
            visualize,
            waveletFunction
            );
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
