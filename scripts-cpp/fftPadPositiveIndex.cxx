#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumberToString.h"
#include "itkViewImage.h"
#include "itkFFTPadPositiveIndexImageFilter.h"
// Boundary Conditions:
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkConstantBoundaryCondition.h"
#include "itkPeriodicBoundaryCondition.h"

// boost::program_options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
namespace po = boost::program_options;
// boost::filesystem
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

template<typename ImageType>
void pad_image(
    const std::string & inputImage,
    const std::string & outputImage,
    const std::string & boundaryConditionName,
    bool visualize)
{
  using ReaderType = itk::ImageFileReader< ImageType >;
  auto reader = ReaderType::New();
  reader->SetFileName( inputImage );
  typename ImageType::Pointer handle = reader->GetOutput();
  using FFTPadFilterType = itk::FFTPadPositiveIndexImageFilter< ImageType >;
  auto fftPadFilter = FFTPadFilterType::New();
  fftPadFilter->SetInput(reader->GetOutput());
  if(boundaryConditionName == "ZeroFluxNeumann")
    {
    // This is the default
    }
  else if(boundaryConditionName == "Constant")
    {
    using BoundaryCondition = itk::ConstantBoundaryCondition< ImageType, ImageType >;
    BoundaryCondition boundaryCondition;
    boundaryCondition.SetConstant(itk::NumericTraits<typename ImageType::PixelType>::ZeroValue());
    fftPadFilter->SetBoundaryCondition(&boundaryCondition);
    }
  else if(boundaryConditionName == "Periodic")
    {
    using BoundaryCondition = itk::PeriodicBoundaryCondition< ImageType, ImageType >;
    BoundaryCondition boundaryCondition;
    fftPadFilter->SetBoundaryCondition(&boundaryCondition);
    }
  else
    {
    std::cerr << "Invalid boundary condition: " << boundaryConditionName << std::endl;
    throw std::runtime_error( "Invalid boundary condition: " + boundaryConditionName);
    }
  fftPadFilter->Update();
  handle = fftPadFilter->GetOutput();

  auto sizeOriginal = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  auto sizeAfterPad = fftPadFilter->GetOutput()->GetLargestPossibleRegion().GetSize();
  std::cout << "Image Padded" << std::endl;
  std::cout << "Original Size:" << sizeOriginal << std::endl;
  std::cout << "After Pad Size:" << sizeAfterPad << std::endl;

  using WriterType = itk::ImageFileWriter< ImageType >;
  auto writer = WriterType::New();
  writer->SetFileName( outputImage );
  writer->SetInput( fftPadFilter->GetOutput() );
  try
    {
    std::cout << "Output in: " << outputImage << std::endl;
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  if(visualize){
    itk::ViewImage<ImageType>::View( reader->GetOutput(), "Input" );
    itk::ViewImage<ImageType>::View( fftPadFilter->GetOutput(), "Padded Output: " + boundaryConditionName );
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
    ( "boundaryCondition,b", po::value<std::string>()->default_value("ZeroFluxNeumann"), "Boundary Condition." )
    ( "visualize,t", po::bool_switch()->default_value(false), "Visualize using vtk based viewer.");

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
  const std::string boundaryConditionName = vm["boundaryCondition"].as<std::string>();
  const bool visualize = vm["visualize"].as<bool>();

  // END PARSE
  itk::NumberToString< double > n2s;
  const fs::path inputImageStem_path = fs::path(inputImage).stem();
  const fs::path outputFolder_path = fs::absolute(fs::path(outputFolder));
  const std::string outputImage = (outputFolder_path /
      fs::path(inputImageStem_path.string() +
        "_fftPad" + boundaryConditionName + "." + outputExtension)).string();

  if(pixelType == "double" ) {
    if(dimension == 3) {
      using ImageType = itk::Image<double, 3>;
      pad_image<ImageType>(inputImage, outputImage, boundaryConditionName, visualize);
    } else if (dimension == 2) {
      using ImageType = itk::Image<double, 2>;
      pad_image<ImageType>(inputImage, outputImage, boundaryConditionName, visualize);
    }
  } else if(pixelType == "float" ) {
    if(dimension == 3) {
      using ImageType = itk::Image<float, 3>;
      pad_image<ImageType>(inputImage, outputImage, boundaryConditionName, visualize);
    } else if (dimension == 2) {
      using ImageType = itk::Image<float, 2>;
      pad_image<ImageType>(inputImage, outputImage, boundaryConditionName, visualize);
    }
  } else if(pixelType == "unsigned char" ) {
    if(dimension == 3) {
      using ImageType = itk::Image<unsigned char, 3>;
      pad_image<ImageType>(inputImage, outputImage, boundaryConditionName, visualize);
    } else if (dimension == 2) {
      using ImageType = itk::Image<unsigned char, 2>;
      pad_image<ImageType>(inputImage, outputImage, boundaryConditionName, visualize);
    }
  } else if(pixelType == "short" ) {
    if(dimension == 3) {
      using ImageType = itk::Image<short, 3>;
      pad_image<ImageType>(inputImage, outputImage, boundaryConditionName, visualize);
    } else if (dimension == 2) {
      using ImageType = itk::Image<short, 2>;
      pad_image<ImageType>(inputImage, outputImage, boundaryConditionName, visualize);
    }
  } else if(pixelType == "int" ) {
    if(dimension == 3) {
      using ImageType = itk::Image<int, 3>;
      pad_image<ImageType>(inputImage, outputImage, boundaryConditionName, visualize);
    } else if (dimension == 2) {
      using ImageType = itk::Image<int, 2>;
      pad_image<ImageType>(inputImage, outputImage, boundaryConditionName, visualize);
    }
  } else {
    std::cout << "PixelType: " << pixelType << " is not supported" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
