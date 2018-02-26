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

// Software Guide : BeginLatex
//
// When the differential equation governing the level set evolution has
// a very simple form, a fast evolution algorithm called fast marching
// can be used.
//
// The following example illustrates the use of the
// \doxygen{FastMarchingImageFilter}. This filter implements a fast marching
// solution to a simple level set evolution problem.  In this example, the
// speed term used in the differential equation is expected to be provided by
// the user in the form of an image.  This image is typically computed as a
// function of the gradient magnitude.  Several mappings are popular in the
// literature, for example, the negative exponential $exp(-x)$ and the
// reciprocal $1/(1+x)$. In the current example we decided to use a Sigmoid
// function since it offers a good number of control parameters that can be
// customized to shape a nice speed image.
//
// The mapping should be done in such a way that the propagation speed of the
// front will be very low close to high image gradients while it will move
// rather fast in low gradient areas. This arrangement will make the contour
// propagate until it reaches the edges of anatomical structures in the image
// and then slow down in front of those edges.  The output of the
// FastMarchingImageFilter is a \emph{time-crossing map} that
// indicates, for each pixel, how much time it would take for the front to
// arrive at the pixel location.
//
// The application of a threshold in the output image is then equivalent to
// taking a snapshot of the contour at a particular time during its evolution.
// It is expected that the contour will take a longer time to cross over
// the edges of a particular anatomical structure. This should result in large
// changes on the time-crossing map values close to the structure edges.
// Segmentation is performed with this filter by locating a time range in which
// the contour was contained for a long time in a region of the image space.
//
// Figure~\ref{fig:FastMarchingCollaborationDiagram} shows the major components
// involved in the application of the FastMarchingImageFilter to a
// segmentation task. It involves an initial stage of smoothing using the
// \doxygen{CurvatureAnisotropicDiffusionImageFilter}. The smoothed image is
// passed as the input to the
// \doxygen{GradientMagnitudeRecursiveGaussianImageFilter} and then to the
// \doxygen{SigmoidImageFilter}.  Finally, the output of the
// FastMarchingImageFilter is passed to a
// \doxygen{BinaryThresholdImageFilter} in order to produce a binary mask
// representing the segmented object.
//
// The code in the following example illustrates the typical setup of a
// pipeline for performing segmentation with fast marching. First, the input
// image is smoothed using an edge-preserving filter. Then the magnitude of its
// gradient is computed and passed to a sigmoid filter. The result of the
// sigmoid filter is the image potential that will be used to affect the speed
// term of the differential equation.
//
// #include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumberToString.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkViewImage.h"

#include "itkRescaleIntensityImageFilter.h"

// boost::program_options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
namespace po = boost::program_options;
// boost::filesystem
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

int main( int argc, char *argv[] )
{
  /*-------------- Parse command line -----------------------------*/
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>()->required(), "Input image." )
    ( "outputFolder,o", po::value<std::string>()->required(), "Outputfolder path. A number of images will be generated depending on levels and bands." )
    ( "outputExtension,e", po::value<std::string>()->default_value("nrrd"), "Output extension." )
    ( "generateGradientMagnitudeOutputImage,x", po::bool_switch()->default_value(false), ".")
    ( "sigma,s", po::value<double>()->required(), "0.5" )
    ( "sigmoidAlpha,a", po::value<double>()->required(), " -0.5." )
    ( "sigmoidBeta,b", po::value<double>()->required(), "2.0." )
    ( "timeThreshold,c", po::value<double>()->required(), "100" )
    ( "stoppingValue,z", po::value<double>()->required(), "100" )
    ( "safeBinaryPercentage,p", po::value<double>()->default_value(0.2), " To generate seeds" )
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
  const double sigma = vm["sigma"].as<double>();
  const double sigmoidAlpha = vm["sigmoidAlpha"].as<double>();
  const double sigmoidBeta = vm["sigmoidBeta"].as<double>();
  const double timeThreshold = vm["timeThreshold"].as<double>();
  const double stoppingValue = vm["stoppingValue"].as<double>();
  const double safeBinaryPercentage = vm["safeBinaryPercentage"].as<double>();
  const bool generateGradientMagnitudeOutputImage =
    vm["generateGradientMagnitudeOutputImage"].as<bool>();
  // END PARSE
  itk::NumberToString< double > n2s;
  const fs::path inputImageStem_path = fs::path(inputImage).stem();
  const fs::path outputFolder_path = fs::absolute(fs::path(outputFolder));
  const std::string parameters = std::string(
      "_s" + n2s(sigma) + "_a" + n2s(sigmoidAlpha) + "_b" + n2s(sigmoidBeta) +
      "_time" + n2s(timeThreshold) + "_stop" + n2s(stoppingValue));
  const std::string outputImage = (outputFolder_path /
      fs::path(inputImageStem_path.string() +
        parameters + "." + outputExtension)).string();

  using InternalPixelType = float;
  const     unsigned int    Dimension = 3;
  using InternalImageType = itk::Image< InternalPixelType, Dimension >;
  // Binary output
  using OutputPixelType = unsigned char;
  using OutputImageType = itk::Image< OutputPixelType, Dimension >;
  using BinaryImageType = OutputImageType;
  using ThresholdingFilterType = itk::BinaryThresholdImageFilter< InternalImageType,
                        OutputImageType    >;
  ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
  //  The upper threshold passed to the BinaryThresholdImageFilter
  //  will define the time snapshot that we are taking from the time-crossing
  //  map. In an ideal application the user should be able to select this
  //  threshold interactively using visual feedback. Here, since it is a
  //  minimal example, the value is taken from the command line arguments.
  thresholder->SetLowerThreshold(           0.0  );
  thresholder->SetUpperThreshold( timeThreshold  );
  thresholder->SetOutsideValue(  0  );
  thresholder->SetInsideValue(  255 );
  using ReaderType = itk::ImageFileReader< InternalImageType >;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImage );
  using WriterType = itk::ImageFileWriter<  OutputImageType  >;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputImage );

  using CastFilterType = itk::RescaleIntensityImageFilter<
                               InternalImageType,
                               OutputImageType >;
  using GradientFilterType = itk::GradientMagnitudeRecursiveGaussianImageFilter<
                               InternalImageType,
                               InternalImageType >;
  GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
  //  The GradientMagnitudeRecursiveGaussianImageFilter performs the
  //  equivalent of a convolution with a Gaussian kernel followed by a
  //  derivative operator. The sigma of this Gaussian can be used to control
  //  the range of influence of the image edges. This filter has been discussed
  //  in Section~\ref{sec:GradientMagnitudeRecursiveGaussianImageFilter}.

  gradientMagnitude->SetSigma(  sigma  );

  using SigmoidFilterType = itk::SigmoidImageFilter<
                               InternalImageType,
                               InternalImageType >;
  SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
  //  The minimum and maximum values of the SigmoidImageFilter output are
  //  defined with the methods \code{SetOutputMinimum()} and
  //  \code{SetOutputMaximum()}. In our case, we want these two values to be
  //  $0.0$ and $1.0$ respectively in order to get a nice speed image to feed
  //  to the FastMarchingImageFilter. Additional details on the use of
  //  the SigmoidImageFilter are presented in
  //  Section~\ref{sec:IntensityNonLinearMapping}.
  sigmoid->SetOutputMinimum(  0.0  );
  sigmoid->SetOutputMaximum(  1.0  );
  //  The SigmoidImageFilter class requires two parameters to define the linear
  //  transformation to be applied to the sigmoid argument. These parameters
  //  are passed using the \code{SetAlpha()} and \code{SetBeta()} methods. In
  //  the context of this example, the parameters are used to intensify the
  //  differences between regions of low and high values in the speed image. In
  //  an ideal case, the speed value should be $1.0$ in the homogeneous regions
  //  of anatomical structures and the value should decay rapidly to $0.0$
  //  around the edges of structures. The heuristic for finding the values is
  //  the following: From the gradient magnitude image, let's call $K1$ the
  //  minimum value along the contour of the anatomical structure to be
  //  segmented. Then, let's call $K2$ an average value of the gradient
  //  magnitude in the middle of the structure. These two values indicate the
  //  dynamic range that we want to map to the interval $[0:1]$ in the speed
  //  image.  We want the sigmoid to map $K1$ to $0.0$ and $K2$ to $1.0$. Given
  //  that $K1$ is expected to be higher than $K2$ and we want to map those
  //  values to $0.0$ and $1.0$ respectively, we want to select a negative
  //  value for alpha so that the sigmoid function will also do an inverse
  //  intensity mapping. This mapping will produce a speed image such that the
  //  level set will march rapidly on the homogeneous region and will
  //  definitely stop on the contour. The suggested value for beta is
  //  $(K1+K2)/2$ while the suggested value for alpha is $(K2-K1)/6$, which
  //  must be a negative number.  In our simple example the values are provided
  //  by the user from the command line arguments. The user can estimate these
  //  values by observing the gradient magnitude image.

  sigmoid->SetAlpha( sigmoidAlpha );
  sigmoid->SetBeta(  sigmoidBeta  );

  using FastMarchingFilterType = itk::FastMarchingImageFilter< InternalImageType,
                              InternalImageType >;
  FastMarchingFilterType::Pointer  fastMarching
                                              = FastMarchingFilterType::New();

  gradientMagnitude->SetInput( reader->GetOutput() );
  sigmoid->SetInput( gradientMagnitude->GetOutput() );
  fastMarching->SetInput( sigmoid->GetOutput() );
  thresholder->SetInput( fastMarching->GetOutput() );
  writer->SetInput( thresholder->GetOutput() );

  //  The FastMarchingImageFilter requires the user to provide a seed point
  //  from which the contour will expand. The user can actually pass not only
  //  one seed point but a set of them. A good set of seed points increases
  //  the chances of segmenting a complex object without missing parts. The
  //  use of multiple seeds also helps to reduce the amount of time needed by
  //  the front to visit a whole object and hence reduces the risk of leaks
  //  on the edges of regions visited earlier. For example, when segmenting
  //  an elongated object, it is undesirable to place a single seed at one
  //  extreme of the object since the front will need a long time to
  //  propagate to the other end of the object. Placing several seeds along
  //  the axis of the object will probably be the best strategy to ensure
  //  that the entire object is captured early in the expansion of the
  //  front. One of the important properties of level sets is their natural
  //  ability to fuse several fronts implicitly without any extra
  //  bookkeeping. The use of multiple seeds takes good advantage of this
  //  property.
  using NodeContainer = FastMarchingFilterType::NodeContainer;
  using NodeType = FastMarchingFilterType::NodeType;
  auto seeds = NodeContainer::New();
  seeds->Initialize();

  using MinMaxCalculator = itk::MinimumMaximumImageCalculator<InternalImageType>;
  auto min_max_calculator = MinMaxCalculator::New();
  reader->Update();
  min_max_calculator->SetImage(reader->GetOutput());
  min_max_calculator->Compute();
  auto max_value = min_max_calculator->GetMaximum();
  auto min_value = min_max_calculator->GetMinimum();
  auto lower_threshold = (max_value - min_value) * safeBinaryPercentage;
  ThresholdingFilterType::Pointer safeBinarizer = ThresholdingFilterType::New();
  //  The upper threshold passed to the BinaryThresholdImageFilter
  //  will define the time snapshot that we are taking from the time-crossing
  //  map. In an ideal application the user should be able to select this
  //  threshold interactively using visual feedback. Here, since it is a
  //  minimal example, the value is taken from the command line arguments.
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
  // safeBinarizer.SetUpperThreshold(max_value);
  // Set seeds from binary image
  using BinaryIteratorType = itk::ImageRegionConstIteratorWithIndex< BinaryImageType >;
  BinaryIteratorType it( safeBinaryOutput, safeBinaryOutput->GetLargestPossibleRegion() );
  // Walk image
  unsigned int count = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
    if (it.Get() > itk::NumericTraits< OutputPixelType >::Zero)
      {
      NodeType node;
      node.SetIndex( it.GetIndex() );
      node.SetValue( 0.0 );
      seeds->InsertElement( count, node );
      count++;
      }
    }

  //  The set of seed nodes is now passed to the FastMarchingImageFilter with
  //  the method \code{SetTrialPoints()}.
  fastMarching->SetTrialPoints(  seeds  );

  //  The FastMarchingImageFilter requires the user to specify the
  //  size of the image to be produced as output. This is done using the
  //  \code{SetOutputSize()} method. Note that the size is obtained here from the
  //  output image of the smoothing filter. The size of this image is valid
  //  only after the \code{Update()} method of this filter has been called
  //  directly or indirectly.
  fastMarching->SetOutputSize(
           reader->GetOutput()->GetBufferedRegion().GetSize() );

  //  Since the front representing the contour will propagate continuously
  //  over time, it is desirable to stop the process once a certain time has
  //  been reached. This allows us to save computation time under the
  //  assumption that the region of interest has already been computed. The
  //  value for stopping the process is defined with the method
  //  \code{SetStoppingValue()}. In principle, the stopping value should be a
  //  little bit higher than the threshold value.
  fastMarching->SetStoppingValue(  stoppingValue  );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  if(generateGradientMagnitudeOutputImage)
  {
    const std::string outputGradientMagnitudeImage =
      (outputFolder_path / fs::path(inputImageStem_path.string() +
                                    parameters + "_GradientMagnitude." +
                                    outputExtension)).string();
    try
    {
      CastFilterType::Pointer caster2 = CastFilterType::New();
      WriterType::Pointer writer2 = WriterType::New();
      caster2->SetInput( gradientMagnitude->GetOutput() );
      writer2->SetInput( caster2->GetOutput() );
      writer2->SetFileName(outputGradientMagnitudeImage);
      caster2->SetOutputMinimum(   0 );
      caster2->SetOutputMaximum( 255 );
      writer2->Update();
    }
    catch( itk::ExceptionObject & err )
    {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  }

  if(visualize){
    itk::Testing::ViewImage( thresholder->GetOutput(), "Output" );
    itk::Testing::ViewImage( gradientMagnitude->GetOutput(), "gradient" );
    itk::Testing::ViewImage( sigmoid->GetOutput(), "sigmoid" );
    itk::Testing::ViewImage( fastMarching->GetOutput(), "FastMarching (no-threshold)" );
  }

  //  \itkcaption[FastMarching segmentation example parameters]{Parameters used
  //  for segmenting some brain structures shown in
  //  Figure~\ref{fig:FastMarchingImageFilterOutput2} using the filter
  //  FastMarchingImageFilter. All of them used a stopping value of
  //  100.\label{tab:FastMarchingImageFilterOutput2}}
  //  \end{table}
  //
  //  Figure~\ref{fig:FastMarchingImageFilterOutput} presents the intermediate
  //  outputs of the pipeline illustrated in
  //  Figure~\ref{fig:FastMarchingCollaborationDiagram}. They are from left to
  //  right: the output of the anisotropic diffusion filter, the gradient
  //  magnitude of the smoothed image and the sigmoid of the gradient magnitude
  //  which is finally used as the speed image for the
  //  FastMarchingImageFilter.
  return EXIT_SUCCESS;
}
