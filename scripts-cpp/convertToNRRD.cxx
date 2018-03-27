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
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkTIFFImageIO.h"
#include "itkImageFileWriter.h"
#include "itkNumberToString.h"
#include "itkViewImage.h"
#include <itkChangeInformationImageFilter.h>

// boost::program_options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
namespace po = boost::program_options;
// boost::filesystem
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

template<typename ImageType>
void convert_image(
    const std::string & inputImage,
    const std::string & outputImage,
    bool changeInfo,
    bool visualize)
{
  using ReaderType = itk::ImageFileReader< ImageType >;
  auto reader = ReaderType::New();
  reader->SetFileName( inputImage );
  // Defense against dodgy tiff images (wrong metadata)
  // Read incorrectly with SCIFIO module.
  auto tiffIO = itk::TIFFImageIO::New();
  if(tiffIO->CanReadFile(inputImage.c_str()))
    reader->SetImageIO( tiffIO );

  typename ImageType::Pointer handle = reader->GetOutput();
  if(changeInfo)
  {
    std::cout << "Changing info" << std::endl;
    using ChangeInformationFilterType = itk::ChangeInformationImageFilter< ImageType >;
    auto changeInputInfoFilter = ChangeInformationFilterType::New();
    typename ImageType::PointType origin_new;
    origin_new.Fill(0);
    typename ImageType::SpacingType spacing_new;
    spacing_new.Fill(1);
    changeInputInfoFilter->SetInput(reader->GetOutput());
    changeInputInfoFilter->ChangeAll();
    changeInputInfoFilter->SetOutputOrigin(origin_new);
    changeInputInfoFilter->SetOutputSpacing(spacing_new);
    changeInputInfoFilter->Update();
    handle = changeInputInfoFilter->GetOutput();
  }
  std::cout << "Changed" << std::endl;

  using WriterType = itk::ImageFileWriter< ImageType >;
  auto writer = WriterType::New();
  writer->SetFileName( outputImage );
  writer->SetInput( handle );
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
    itk::Testing::ViewImage( reader->GetOutput(), "Input" );
  }

}


int main( int argc, char *argv[])
{
  /*-------------- Parse command line -----------------------------*/
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>()->required(), "Input image." )
    ( "outputFolder,o", po::value<std::string>()->required(), "Outputfolder path. A number of images will be generated depending on levels and bands." )
    ( "outputExtension,e", po::value<std::string>()->default_value("nrrd"), "Output extension." )
    ( "dimension,d", po::value<unsigned int>()->required(), "Use imageInfo if needed." )
    ( "pixelType,p", po::value<std::string>()->required(), "Use imageInfo if needed." )
    ( "visualize,t", po::bool_switch()->default_value(false), "Visualize using vtk based viewer.")
    ( "changeInfo,c",  po::bool_switch()->default_value(false), "ChangeInfo to default." );

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
  const unsigned int dimension = vm["dimension"].as<unsigned int>();
  const std::string pixelType = vm["pixelType"].as<std::string>();
  const bool changeInfo = vm["changeInfo"].as<bool>();
  const bool visualize = vm["visualize"].as<bool>();

  // END PARSE
  itk::NumberToString< double > n2s;
  const fs::path inputImageStem_path = fs::path(inputImage).stem();
  const fs::path outputFolder_path = fs::absolute(fs::path(outputFolder));
  const std::string outputImage = (outputFolder_path /
      fs::path(inputImageStem_path.string() + "." + outputExtension)).string();

  if(pixelType == "double" ) {
    if(dimension == 3) {
      using ImageType = itk::Image<double, 3>;
      convert_image<ImageType>(inputImage, outputImage, changeInfo, visualize);
    } else if (dimension == 2) {
      using ImageType = itk::Image<double, 2>;
      convert_image<ImageType>(inputImage, outputImage, changeInfo, visualize);
    }
  } else if(pixelType == "float" ) {
    if(dimension == 3) {
      using ImageType = itk::Image<float, 3>;
      convert_image<ImageType>(inputImage, outputImage, changeInfo, visualize);
    } else if (dimension == 2) {
      using ImageType = itk::Image<float, 2>;
      convert_image<ImageType>(inputImage, outputImage, changeInfo, visualize);
    }
  } else if(pixelType == "unsigned char" ) {
    if(dimension == 3) {
      using ImageType = itk::Image<unsigned char, 3>;
      convert_image<ImageType>(inputImage, outputImage, changeInfo, visualize);
    } else if (dimension == 2) {
      using ImageType = itk::Image<unsigned char, 2>;
      convert_image<ImageType>(inputImage, outputImage, changeInfo, visualize);
    }
  } else if(pixelType == "unsigned short" ) {
    if(dimension == 3) {
      using ImageType = itk::Image<unsigned short, 3>;
      convert_image<ImageType>(inputImage, outputImage, changeInfo, visualize);
    } else if (dimension == 2) {
      using ImageType = itk::Image<unsigned short, 2>;
      convert_image<ImageType>(inputImage, outputImage, changeInfo, visualize);
    }
  }

  return EXIT_SUCCESS;
}
