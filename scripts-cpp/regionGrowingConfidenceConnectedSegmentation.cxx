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

//  Software Guide : BeginCommandLineArgs
//    INPUTS:  {BrainProtonDensitySlice.png}
//    OUTPUTS: {ConnectedThresholdOutput1.png}
//    ARGUMENTS:    60 116 150 180
//  Software Guide : EndCommandLineArgs
//  Software Guide : BeginCommandLineArgs
//    INPUTS:  {BrainProtonDensitySlice.png}
//    OUTPUTS: {ConnectedThresholdOutput2.png}
//    ARGUMENTS:    81 112 210 250
//  Software Guide : EndCommandLineArgs
//  Software Guide : BeginCommandLineArgs
//    INPUTS:  {BrainProtonDensitySlice.png}
//    OUTPUTS: {ConnectedThresholdOutput3.png}
//    ARGUMENTS:    107 69 180 210
#include "itkConfidenceConnectedImageFilter.h"
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkImageFileReader.h"
#include "itkTIFFImageIO.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumberToString.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkViewImage.h"

// boost::program_options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
namespace po = boost::program_options;
// boost::filesystem
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

int main( int argc, char *argv[])
{
  /*-------------- Parse command line -----------------------------*/
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>()->required(), "Input image." )
    ( "outputFolder,o", po::value<std::string>()->required(), "Outputfolder path. A number of images will be generated depending on levels and bands." )
    ( "outputExtension,e", po::value<std::string>()->default_value("nrrd"), "Output extension." )
    ( "iters,n", po::value<unsigned int>()->required(), "5" )
    ( "multiplier,m", po::value<double>()->required(), "2.5" )
    ( "radius,r", po::value<unsigned int>()->default_value(2), "2" )
    ( "safeBinaryPercentage,p", po::value<double>()->default_value(0.1), " To generate seeds" )

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
  const std::string outputFolder = vm["outputFolder"].as<std::string>();
  const std::string outputExtension = vm["outputExtension"].as<std::string>();
  const bool verbose = vm["verbose"].as<bool>();
  const bool visualize = vm["visualize"].as<bool>();
  const double safeBinaryPercentage = vm["safeBinaryPercentage"].as<double>();
  const unsigned int iters = vm["iters"].as<unsigned int>();
  const double multiplier = vm["multiplier"].as<double>();
  const unsigned int radius = vm["radius"].as<unsigned int>();

  // END PARSE
  itk::NumberToString< double > n2s;
  const fs::path inputImageStem_path = fs::path(inputImage).stem();
  const fs::path outputFolder_path = fs::absolute(fs::path(outputFolder));
  const std::string parameters =
    "_SegmentRG_CC_iters" + n2s(iters) +
    "_m" + n2s(multiplier) +
    "_r" + n2s(radius) +
    "_p" + n2s(safeBinaryPercentage);
  const std::string outputImage = (outputFolder_path /
      fs::path(inputImageStem_path.string() +
        parameters + "." + outputExtension)).string();

  using InternalPixelType = float;
  const     unsigned int    Dimension = 3;
  using InternalImageType = itk::Image< InternalPixelType, Dimension >;

  using OutputPixelType = unsigned char;
  using OutputImageType = itk::Image< OutputPixelType, Dimension >;
  using CastingFilterType =
      itk::CastImageFilter< InternalImageType, OutputImageType >;
  auto caster = CastingFilterType::New();

  using ReaderType = itk::ImageFileReader< InternalImageType >;
  using WriterType = itk::ImageFileWriter<  OutputImageType  >;

  auto reader = ReaderType::New();
  reader->SetFileName( inputImage );
  // Defense against dodgy tiff images (wrong metadata)
  // Read incorrectly with SCIFIO module.
  auto tiffIO = itk::TIFFImageIO::New();
  if(tiffIO->CanReadFile(inputImage.c_str()))
    reader->SetImageIO( tiffIO );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  auto writer = WriterType::New();
  writer->SetFileName( outputImage );

  using ConnectedFilterType =
    itk::ConfidenceConnectedImageFilter< InternalImageType, InternalImageType >;
  auto connectedFilter = ConnectedFilterType::New();
  connectedFilter->SetInput( reader->GetOutput() );
  caster->SetInput( connectedFilter->GetOutput() );
  writer->SetInput( caster->GetOutput() );

  connectedFilter->SetReplaceValue( 255 );
  connectedFilter->SetMultiplier( multiplier );
  connectedFilter->SetNumberOfIterations( iters );
  connectedFilter->SetInitialNeighborhoodRadius( radius );

  // Add seeds based on binarization:
  using MinMaxCalculator = itk::MinimumMaximumImageCalculator<InternalImageType>;
  auto min_max_calculator = MinMaxCalculator::New();
  min_max_calculator->SetImage(reader->GetOutput());
  min_max_calculator->Compute();
  auto max_value = min_max_calculator->GetMaximum();
  auto min_value = min_max_calculator->GetMinimum();
  auto lower_threshold = (max_value - min_value) * safeBinaryPercentage;
  using ThresholdingFilterType =
    itk::BinaryThresholdImageFilter< InternalImageType, OutputImageType >;
  auto safeBinarizer = ThresholdingFilterType::New();
  safeBinarizer->SetInput(reader->GetOutput());
  // safeBinarizer->SetOutsideValue(  0  );
  // safeBinarizer->SetInsideValue(  255 );
  safeBinarizer->SetLowerThreshold(lower_threshold);
  safeBinarizer->SetUpperThreshold(max_value);
  safeBinarizer->Update();
  if(verbose){
    std::cout << "Min: " << min_value << std::endl;
    std::cout << "Max: " << max_value << std::endl;
    std::cout << "lower_threshold " << lower_threshold << std::endl;
  }
  auto safeBinaryOutput = safeBinarizer->GetOutput();

  if(visualize){
    itk::Testing::ViewImage( reader->GetOutput(), "Input" );
    itk::Testing::ViewImage( safeBinaryOutput, "Safe Binary Output" );
  }
  // Set seeds from binary image
  using OutputIteratorType = itk::ImageRegionConstIteratorWithIndex< OutputImageType >;
  OutputIteratorType it( safeBinaryOutput, safeBinaryOutput->GetLargestPossibleRegion() );
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
    if (it.Get() > itk::NumericTraits< OutputPixelType >::Zero)
      {
      connectedFilter->AddSeed( it.GetIndex() );
      }
    }

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  if(visualize){
    itk::Testing::ViewImage( caster->GetOutput(), "Output" );
  }


  return EXIT_SUCCESS;
}
