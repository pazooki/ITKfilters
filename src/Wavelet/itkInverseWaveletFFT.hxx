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
#ifndef itkInverseWaveletFT_hxx
#define itkInverseWaveletFT_hxx
#include "visualize_functions.h" // TODO REMOVE
#include <array>
#include <itkSubtractImageFilter.h>
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include "itkInverseWaveletFT.h"
#include <itkComposeImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkVectorImageToImageAdaptor.h>
#include "itkMultiplyImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkProgressAccumulator.h"
#include <itkAddImageFilter.h>
#include <itkSquareImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>//For test max/min before sum.
#include <itkSqrtImageFilter.h>
#include <itkAtanImageFilter.h> // For phase studies
#include <itkDivideImageFilter.h>
// Eigen Calculations
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkNeighborhoodInnerProduct.h>
#include <itkGaussianImageSource.h>
#include <itkMatrix.h>
#include <itkSymmetricEigenAnalysis.h>

namespace itk
{
template< typename TInputImage >
InverseWaveletFT< TInputImage >
::InverseWaveletFT()
{
  this->SetNumberOfRequiredOutputs(4);
  this->SetNumberOfRequiredInputs(1);

  this->SetNthOutput( 0, this->MakeOutput(0) );
  this->SetNthOutput( 1, this->MakeOutput(1) );
  this->SetNthOutput( 2, this->MakeOutput(2) );
  this->SetNthOutput( 3, this->MakeOutput(3) );

  this->m_SigmaGaussianDerivative = 1.0;
  this->m_StatisticsMean = 0.0;
}

template< typename TInputImage >
void
InverseWaveletFT< TInputImage >
::GenerateData()
{

  typename StatisticsImageFilterType::Pointer statisticsImageFilter
          = StatisticsImageFilterType::New ();
  statisticsImageFilter->SetInput(this->GetInput());
  statisticsImageFilter->Update();
  this->m_StatisticsMean = statisticsImageFilter->GetMean();
  // // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);
  // Allocate the outputs
  InputImageType* evenPart = this->GetOutputReal();
  evenPart->SetRegions( evenPart->GetLargestPossibleRegion() );
  evenPart->Allocate();
  RieszComponentsImageType* rieszComponents = this->GetOutputRieszComponents();
  rieszComponents->SetRegions( rieszComponents->GetRequestedRegion() );
  rieszComponents->Allocate();
  InputImageType* normPart = this->GetOutputRieszNorm();
  normPart->SetRegions( normPart->GetRequestedRegion() );
  normPart->Allocate();

  //Substract Mean value of image to every pixel to remove DC comp
  typedef itk::SubtractImageFilter<InputImageType> SubtractFilterType;
  typename SubtractFilterType::Pointer subtractFilter =
    SubtractFilterType::New();
  subtractFilter->SetInput(this->GetInput());
  subtractFilter->SetConstant(this->m_StatisticsMean);
  subtractFilter->Update();

  typename FFTFilterType::Pointer fftFilter = FFTFilterType::New();
  // fftFilter->SetInput( this->GetInput() );
  fftFilter->SetInput( subtractFilter->GetOutput());
  progress->RegisterInternalFilter(fftFilter, 1.0f);
  fftFilter->Update();
  std::cout << fftFilter->GetOutput()->GetBufferedRegion() << std::endl;
  typename ComplexImageType::Pointer fftForward = this->GetOutputFFT();
  fftForward->SetRegions(fftFilter->GetOutput()->GetRequestedRegion());
  fftForward->Allocate();
  // fftForward = fftFilter->GetOutput();
  ImageAlgorithm::Copy(fftFilter->GetOutput(), fftForward.GetPointer(),
      fftFilter->GetOutput()->GetRequestedRegion(), fftFilter->GetOutput()->GetRequestedRegion());

  // this->PrintSelf(std::cout,2);

  std::cout << "RealComponent" << std::endl;
  InputImagePointer realC = this->ComputeRealComponent(fftForward);
  this->GetOutputReal()->Graft(realC);
  std::cout << "RieszComponents" << std::endl;
  // RIESZ OUTPUTS
  typename RieszComponentsImageType::Pointer rieszC =
    this->ComputeRieszComponents(fftForward);
  this->GetOutputRieszComponents()->Graft(rieszC);
  std::cout << "RieszNorm" << std::endl;
  InputImagePointer norm = this->ComputeRieszNorm(rieszC);
  this->GetOutputRieszNorm()->Graft(norm);
  std::cout << "End GenerateData" << std::endl;
}

template< typename TInputImage >
void
InverseWaveletFT< TInputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Sigma of Gaussian Derivative : " << m_SigmaGaussianDerivative
     << std::endl;
}

template< typename TInputImage>
DataObject::Pointer InverseWaveletFT<TInputImage>
::MakeOutput(DataObjectPointerArraySizeType idx)
{
  DataObject::Pointer output;

  switch ( idx )
    {
    case 0:
      output = ( InputImageType::New() ).GetPointer();
      break;
    case 1:
      output = ( RieszComponentsImageType::New() ).GetPointer();
      break;
    case 2:
      output = ( InputImageType::New() ).GetPointer();
      break;
    case 3:
      output = ( ComplexImageType::New() ).GetPointer();
      break;
    default:
      std::cerr << "No output " << idx << std::endl;
      output = NULL;
      break;
    }
  return output.GetPointer();
}

} // end namespace itk
#endif
