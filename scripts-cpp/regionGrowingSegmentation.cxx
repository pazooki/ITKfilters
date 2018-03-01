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
#include "itkConnectedThresholdImageFilter.h"
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkImageFileReader.h"
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

  // if( argc < 7 )
  //   {
  //   std::cerr << "Missing Parameters " << std::endl;
  //   std::cerr << "Usage: " << argv[0];
  //   std::cerr << " inputImage  outputImage seedX seedY lowerThreshold upperThreshold" << std::endl;
  //   return EXIT_FAILURE;
  //   }
int main( int argc, char *argv[])
{
  /*-------------- Parse command line -----------------------------*/
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>()->required(), "Input image." )
    ( "outputFolder,o", po::value<std::string>()->required(), "Outputfolder path. A number of images will be generated depending on levels and bands." )
    ( "outputExtension,e", po::value<std::string>()->default_value("nrrd"), "Output extension." )
    ( "lowerThreshold,l", po::value<double>()->required(), "0.5" )
    ( "upperThreshold,u", po::value<double>()->required(), " -0.5." )
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
  const double lowerThreshold = vm["lowerThreshold"].as<double>();
  const double upperThreshold = vm["upperThreshold"].as<double>();
  // END PARSE
  itk::NumberToString< double > n2s;
  const fs::path inputImageStem_path = fs::path(inputImage).stem();
  const fs::path outputFolder_path = fs::absolute(fs::path(outputFolder));
  const std::string parameters =
      "_SegmentRG_l" + n2s(lowerThreshold) +
      "_u" + n2s(upperThreshold) +
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

  using ConnectedFilterType = itk::ConnectedThresholdImageFilter< InternalImageType,
                                    InternalImageType >;
  auto connectedThreshold = ConnectedFilterType::New();
  connectedThreshold->SetInput( reader->GetOutput() );
  caster->SetInput( connectedThreshold->GetOutput() );
  writer->SetInput( caster->GetOutput() );

  connectedThreshold->SetLower(  lowerThreshold  );
  connectedThreshold->SetUpper(  upperThreshold  );
  connectedThreshold->SetReplaceValue( 255 );

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
      connectedThreshold->AddSeed( it.GetIndex() );
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
