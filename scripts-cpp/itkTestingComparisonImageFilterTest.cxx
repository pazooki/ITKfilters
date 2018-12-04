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

#include <iostream>
#include "itkTestingComparisonImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <cstdlib>

template< typename InputImageType, typename OutputImageType>
int compareImages(
    const std::string& input1Valid,
    const std::string& input2Test,
    const std::string& outputDifference,
    int differenceThreshold,
    int toleranceRadius,
    unsigned long allowedNumberOfPixelsWithDifferences
    )
{
  bool inputsAreSimilar = true;
  using ReaderType = itk::ImageFileReader< InputImageType  >;

  auto reader1 = ReaderType::New();
  auto reader2 = ReaderType::New();

  reader1->SetFileName( input1Valid );
  reader2->SetFileName( input2Test );

  // Define the filter
  using FilterType = itk::Testing::ComparisonImageFilter<
                             InputImageType,
                             OutputImageType >;

  auto filter = FilterType::New();

  // setup the filter
  filter->SetDifferenceThreshold( differenceThreshold );
  filter->SetToleranceRadius( toleranceRadius );

  // wire the pipeline
  filter->SetValidInput( reader1->GetOutput() );
  filter->SetTestInput(  reader2->GetOutput() );

  // Write the output
  using WriterType = itk::ImageFileWriter< OutputImageType >;

  auto writer = WriterType::New();

  writer->SetInput( filter->GetOutput() );

  writer->SetFileName( outputDifference );

  writer->Update();

  unsigned long currentNumberOfPixelsWithDifferences =
    filter->GetNumberOfPixelsWithDifferences();

  if(currentNumberOfPixelsWithDifferences  > allowedNumberOfPixelsWithDifferences)
  {
    inputsAreSimilar = false;
    std::cout << "FAIL. Images not similar.\n  CurrentNumberOfPixelsWithDifferences ( " <<
      currentNumberOfPixelsWithDifferences << " ) greater than allowed ( " <<
      allowedNumberOfPixelsWithDifferences << " )." << std::endl;
  }
  if (inputsAreSimilar)
  {
    std::cout << "SUCCESS. Images are similar." << std::endl;
  }

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  if( argc < 6  ||  argc > 10)
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0];
    std::cerr <<
      "  inputImageFile1 inputImageFile2 outputImage DimensionOfImages inputPixelType\n"
      " [outputPixelType = inputPixelType] [threshold = 0] [radius = 0] [numberOfPixelsWithDifferences = 0]\n"
      "Hint: You might want to use imageInfo to get the dimension and pixel types."
      << std::endl;
    return EXIT_FAILURE;
    }

  std::string input1Valid(argv[1]);
  std::string input2Test(argv[2]);
  std::string outputDifference(argv[3]);
  unsigned int dimension = std::stoi(argv[4]);
  std::string inputPixelType(argv[5]);
  // optionals
  std::string outputPixelType = inputPixelType;
  int differenceThreshold = 0;
  int toleranceRadius = 0;
  unsigned long allowedNumberOfPixelsWithDifferences = 0;
  if( argc > 6 )
    outputPixelType = argv[6];
  if( argc > 7 )
    differenceThreshold = std::stoi( argv[7] );
  if( argc > 8 )
    toleranceRadius = std::stoi( argv[8] );
  if( argc > 9 )
    unsigned long allowedNumberOfPixelsWithDifferences = std::stoul(argv[9]);

  std::cout << "--- Input Parameters ---" << input1Valid << std::endl;
  std::cout << "Input1: Valid Image: " << input1Valid << std::endl;
  std::cout << "Input2: Test Image: "  << input2Test << std::endl;
  std::cout << "Output: Difference Image: "  << outputDifference << std::endl;
  std::cout << "DifferenceThreshold:" << differenceThreshold << std::endl;
  std::cout << "ToleranceRadius: "  << toleranceRadius << std::endl;
  std::cout << "Allowed number of pixels with differences: "  <<
    allowedNumberOfPixelsWithDifferences << std::endl;

  if(inputPixelType == "double" ) {
    if(dimension == 3) {
      using InputImageType = itk::Image<double, 3>;
      using OutputImageType = itk::Image<double, 3>;
      compareImages<InputImageType, OutputImageType>(
          input1Valid, input2Test, outputDifference, differenceThreshold,
          toleranceRadius, allowedNumberOfPixelsWithDifferences);
    } else if (dimension == 2) {
      using InputImageType = itk::Image<double, 2>;
      using OutputImageType = itk::Image<double, 2>;
      compareImages<InputImageType, OutputImageType>(
          input1Valid, input2Test, outputDifference, differenceThreshold,
          toleranceRadius, allowedNumberOfPixelsWithDifferences);
    }
  } else if(inputPixelType == "float" ) {
    if(dimension == 3) {
      using InputImageType = itk::Image<float, 3>;
      using OutputImageType = itk::Image<float, 3>;
      compareImages<InputImageType, OutputImageType>(
          input1Valid, input2Test, outputDifference, differenceThreshold,
          toleranceRadius, allowedNumberOfPixelsWithDifferences);
    } else if (dimension == 2) {
      using InputImageType = itk::Image<float, 2>;
      using OutputImageType = itk::Image<float, 2>;
      compareImages<InputImageType, OutputImageType>(
          input1Valid, input2Test, outputDifference, differenceThreshold,
          toleranceRadius, allowedNumberOfPixelsWithDifferences);
    }
  } else if(inputPixelType == "unsigned char" ) {
    if(dimension == 3) {
      using InputImageType = itk::Image<unsigned char, 3>;
      using OutputImageType = itk::Image<unsigned char, 3>;
      compareImages<InputImageType, OutputImageType>(
          input1Valid, input2Test, outputDifference, differenceThreshold,
          toleranceRadius, allowedNumberOfPixelsWithDifferences);
    } else if (dimension == 2) {
      using InputImageType = itk::Image<unsigned char, 2>;
      using OutputImageType = itk::Image<unsigned char, 2>;
      compareImages<InputImageType, OutputImageType>(
          input1Valid, input2Test, outputDifference, differenceThreshold,
          toleranceRadius, allowedNumberOfPixelsWithDifferences);
    }
  } else if(inputPixelType == "short" ) {
    if(dimension == 3) {
      using InputImageType = itk::Image<short, 3>;
      using OutputImageType = itk::Image<short, 3>;
      compareImages<InputImageType, OutputImageType>(
          input1Valid, input2Test, outputDifference, differenceThreshold,
          toleranceRadius, allowedNumberOfPixelsWithDifferences);
    } else if (dimension == 2) {
      using InputImageType = itk::Image<short, 2>;
      using OutputImageType = itk::Image<short, 2>;
      compareImages<InputImageType, OutputImageType>(
          input1Valid, input2Test, outputDifference, differenceThreshold,
          toleranceRadius, allowedNumberOfPixelsWithDifferences);
    }
  } else if(inputPixelType == "int" ) {
    if(dimension == 3) {
      using InputImageType = itk::Image<int, 3>;
      using OutputImageType = itk::Image<int, 3>;
      compareImages<InputImageType, OutputImageType>(
          input1Valid, input2Test, outputDifference, differenceThreshold,
          toleranceRadius, allowedNumberOfPixelsWithDifferences);
    } else if (dimension == 2) {
      using InputImageType = itk::Image<int, 2>;
      using OutputImageType = itk::Image<int, 2>;
      compareImages<InputImageType, OutputImageType>(
          input1Valid, input2Test, outputDifference, differenceThreshold,
          toleranceRadius, allowedNumberOfPixelsWithDifferences);
    }
  } else {
    std::cout << "inpuPixelType: " << inputPixelType << " is not supported" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
