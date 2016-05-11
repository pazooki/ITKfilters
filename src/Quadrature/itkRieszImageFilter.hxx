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
#ifndef itkRieszImageFilter_hxx
#define itkRieszImageFilter_hxx
#include "visualize_functions.h" // TODO REMOVE
#include <array>
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include "itkRieszImageFilter.h"
#include <itkComposeImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkVectorImageToImageAdaptor.h>
#include "itkRealToHalfHermitianForwardFFTImageFilter.h"
#include "itkHalfHermitianToRealInverseFFTImageFilter.h"
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
RieszImageFilter< TInputImage >
::RieszImageFilter()
{
  this->SetNumberOfRequiredOutputs(3);
  this->SetNumberOfRequiredInputs(1);

  this->SetNthOutput( 0, this->MakeOutput(0) );
  this->SetNthOutput( 1, this->MakeOutput(1) );
  this->SetNthOutput( 2, this->MakeOutput(2) );

  this->m_SigmaGaussianDerivative = 1.0;
}

template< typename TInputImage >
void
RieszImageFilter< TInputImage >
::GenerateData()
{

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

  typedef RealToHalfHermitianForwardFFTImageFilter < InputImageType> FFTFilterType;
  typename FFTFilterType::Pointer fftFilter = FFTFilterType::New();
  fftFilter->SetInput( this->GetInput() );
  progress->RegisterInternalFilter(fftFilter, 1.0f);
  fftFilter->Update();
  std::cout << fftFilter->GetOutput()->GetBufferedRegion() << std::endl;

  // Save the size and spacing of fftForward in non-accesible members.
  this->m_inputSizeSquare    = fftFilter->GetOutput()->GetBufferedRegion().GetSize() * fftFilter->GetOutput()->GetBufferedRegion().GetSize();
  // Spacing is a vector, and * perform inner product.
  this->m_inputSpacingSquare = fftFilter->GetOutput()->GetSpacing();
  for (unsigned int i = 0; i < ImageDimension ; i++)
    this->m_inputSpacingSquare[i] *= this->m_inputSpacingSquare[i];
  // this->PrintSelf(std::cout,2);

  std::cout << "RealComponent" << std::endl;
  InputImagePointer realC = this->ComputeRealComponent(fftFilter->GetOutput());
  this->GetOutputReal()->Graft(realC);
  std::cout << "RieszComponents" << std::endl;
  // RIESZ OUTPUTS
  typename RieszComponentsImageType::Pointer rieszC = this->ComputeRieszComponents(fftFilter->GetOutput());
  this->GetOutputRieszComponents()->Graft(rieszC);
  std::cout << "RieszNorm" << std::endl;
  InputImagePointer norm = this->ComputeRieszNorm(rieszC);
  this->GetOutputRieszNorm()->Graft(norm);
  std::cout << "End GenerateData" << std::endl;
}

template< typename TInputImage>
typename TInputImage::Pointer
RieszImageFilter<TInputImage>
::ComputeRealComponent(const ComplexImageType* fftForward ) const
{
  typedef itk::ImageDuplicator< ComplexImageType > DuplicatorType;
  typename DuplicatorType::Pointer duplicator= DuplicatorType::New();
  duplicator->SetInputImage(fftForward);
  duplicator->Update();
  typename ComplexImageType::Pointer intermediatePtr = duplicator->GetModifiableOutput();
  intermediatePtr->DisconnectPipeline();
  ImageRegionIteratorWithIndex< ComplexImageType >
  outIt(intermediatePtr, intermediatePtr->GetRequestedRegion() );
  for ( outIt.GoToBegin() ; !outIt.IsAtEnd(); ++outIt )
    {
    const typename ComplexImageType::IndexType index = outIt.GetIndex();
    typename ComplexImageType::PointType evalPoint; // w = (u,v,..);
    intermediatePtr->TransformIndexToPhysicalPoint(index, evalPoint);
    RealType w2 = 0;
    for (unsigned int i = 0; i < ImageDimension ; i++)
      w2 += ( m_inputSpacingSquare[i] / static_cast<RealType>(m_inputSizeSquare[i])) * evalPoint[i] * evalPoint[i];
    const InputImagePixelType w_mod = sqrt( w2 );

    const ComplexImagePixelType value = outIt.Get() *
      w_mod * std::exp(- this->m_SigmaGaussianDerivative * this->m_SigmaGaussianDerivative * w2 ) ;

    outIt.Set( value );
    }
  intermediatePtr->Modified();

  typedef HalfHermitianToRealInverseFFTImageFilter < ComplexImageType, InputImageType> InverseFFTFilterType;
  typename InverseFFTFilterType::Pointer inverseFilter = InverseFFTFilterType::New();
  inverseFilter->SetInput( intermediatePtr );
  inverseFilter->Update();
  typename InputImageType::Pointer output = inverseFilter->GetOutput();
  output->DisconnectPipeline();
  std::cout << output->GetBufferedRegion() << std::endl;
  return output;

}

template< typename TInputImage>
typename RieszImageFilter<TInputImage>::InputImageType::Pointer
RieszImageFilter<TInputImage>
::ComputeRieszComponent(const ComplexImageType* fftForward, const unsigned int & NComponent) const
{

  typedef itk::ImageDuplicator< ComplexImageType > DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(fftForward);
  duplicator->Update();
  typename ComplexImageType::Pointer componentImage = duplicator->GetOutput();
  componentImage->DisconnectPipeline();

  ImageRegionIteratorWithIndex< ComplexImageType >
    outIt( componentImage, componentImage->GetBufferedRegion() );
  for ( outIt.GoToBegin() ; !outIt.IsAtEnd(); ++outIt )
  {
    const typename ComplexImageType::IndexType index = outIt.GetIndex();
    typename ComplexImageType::PointType evalPoint; // w = (u,v,..);
    componentImage->TransformIndexToPhysicalPoint(index, evalPoint);
    RealType w2 = 0;
    for (unsigned int i = 0; i < ImageDimension ; i++)
      w2 += ( m_inputSpacingSquare[i] / static_cast<RealType>(m_inputSizeSquare[i])) * evalPoint[i] * evalPoint[i];
    const InputImagePixelType w_mod = sqrt( w2 );
    const InputImagePixelType shared_value = std::exp( - this->m_SigmaGaussianDerivative * this->m_SigmaGaussianDerivative * w2) ;
    const ComplexImagePixelType value = outIt.Get() *
      std::complex<InputImagePixelType>(w_mod * shared_value,
          -shared_value * static_cast<InputImagePixelType>(evalPoint[NComponent]));
    // std::cout <<"Index: "<< index << " ; outIt.Get(): " << outIt.Get() << std::endl;
    // std::cout <<"wmod: " << w_mod  <<" ; out1: real: " << value.real() << " ccomplex:" << value.imag() <<  std::endl;
    outIt.Set( value );
  }
  componentImage->Modified();

  typedef HalfHermitianToRealInverseFFTImageFilter < ComplexImageType, InputImageType> InverseFFTFilterType;
  typename InverseFFTFilterType::Pointer inverseFilter = InverseFFTFilterType::New();
  inverseFilter->SetInput( componentImage );
  inverseFilter->Update();
  typename InputImageType::Pointer output = inverseFilter->GetOutput();
  output->DisconnectPipeline();
  return output;
}

template< typename TInputImage>
typename RieszImageFilter<TInputImage>::RieszComponentsImageType::Pointer
RieszImageFilter<TInputImage>
::ComputeRieszComponents(const ComplexImageType* fftForward) const
{
  typedef itk::ComposeImageFilter<InputImageType> ComposerType;
  typename ComposerType::Pointer composer= ComposerType::New();
  for(unsigned int N = 0 ; N < ImageDimension ; ++N)
    composer->SetInput(N, ComputeRieszComponent(fftForward, N));
  composer->Update();

  typename RieszComponentsImageType::Pointer output = composer->GetOutput();
  output->DisconnectPipeline();
  return output;
}

template< typename TInputImage>
typename RieszImageFilter<TInputImage>::InputImageType::Pointer
RieszImageFilter<TInputImage>
::ComputeLocalAmplitude(const InputImageType* real_part, const InputImageType* riesz_norm_part) const
{
  typedef itk::SquareImageFilter<InputImageType,InputImageType> SquareImageFilterType;
  auto squareFilterEven = SquareImageFilterType::New();
  squareFilterEven->SetInput(real_part);
  squareFilterEven->Update();
  auto squareFilterNorm = SquareImageFilterType::New();
  squareFilterNorm->SetInput(riesz_norm_part);
  squareFilterNorm->Update();

  typedef itk::AddImageFilter<InputImageType,InputImageType, InputImageType> AddImageFilterType;
  auto addFilter = AddImageFilterType::New();
  addFilter->SetInput1(squareFilterEven->GetOutput());
  addFilter->SetInput2(squareFilterNorm->GetOutput());
  addFilter->Update();

  typedef SqrtImageFilter<InputImageType,InputImageType> SqrtImageFilterType;
  typename SqrtImageFilterType::Pointer sqrtFilter = SqrtImageFilterType::New();
  sqrtFilter->SetInput(addFilter->GetOutput());
  sqrtFilter->Update();
  return sqrtFilter->GetOutput();
}

template< typename TInputImage>
typename RieszImageFilter<TInputImage>::InputImageType::Pointer
RieszImageFilter<TInputImage>
::ComputeRieszWeightedNorm(const RieszComponentsImageType* rieszComponents, const DirectionType & weights) const
{
  typename InputImageType::Pointer output = InputImageType::New();
  output->SetRegions(rieszComponents->GetBufferedRegion());
  output->Allocate();
  output->FillBuffer(NumericTraits<InputImagePixelType>::Zero);

  typedef VectorIndexSelectionCastImageFilter<RieszComponentsImageType,InputImageType> CastIndexType;
  typename CastIndexType::Pointer castIndex = CastIndexType::New();
  castIndex->SetInput(rieszComponents);

  typedef SquareImageFilter<InputImageType,InputImageType> SquareImageFilterType;
  typename SquareImageFilterType::Pointer squareFilter = SquareImageFilterType::New();

  typedef AddImageFilter<InputImageType,InputImageType,InputImageType> AddImageFilterType;
  typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();

  typedef itk::MultiplyImageFilter<TInputImage,TInputImage,TInputImage> MultiplyImageFilterType;
  typename MultiplyImageFilterType::Pointer multiplyByConstant = MultiplyImageFilterType::New();
  for (unsigned int i = 0 ; i<ImageDimension ; ++i){
    castIndex->SetIndex(i);
    castIndex->Update();
    multiplyByConstant->SetInput(castIndex->GetOutput());
    multiplyByConstant->SetConstant(weights[i]);
    multiplyByConstant->Update();
    squareFilter->SetInput(multiplyByConstant->GetOutput());
    squareFilter->Update();
    typename InputImageType::Pointer squareResult = squareFilter->GetOutput();
    squareResult->DisconnectPipeline();
    addFilter->SetInput1(squareResult);
    addFilter->SetInput2(output);
    addFilter->Update();
    output = addFilter->GetOutput();
    output->DisconnectPipeline();
  }
  typedef SqrtImageFilter<InputImageType,InputImageType> SqrtImageFilterType;
  typename SqrtImageFilterType::Pointer sqrtFilter = SqrtImageFilterType::New();
  sqrtFilter->SetInput(output);
  sqrtFilter->InPlaceOn();
  sqrtFilter->Update();
  return sqrtFilter->GetOutput();
}

template< typename TInputImage>
typename RieszImageFilter<TInputImage>::InputImageType::Pointer
RieszImageFilter<TInputImage>
::ComputeRieszNorm(const RieszComponentsImageType* rieszComponents) const
{
  DirectionType weights;
  for (unsigned int i = 0 ; i < ImageDimension ; ++i)
    weights[i] = 1.0 ;
  return ComputeRieszWeightedNorm(rieszComponents, weights);
}

template< typename TInputImage>
typename RieszImageFilter<TInputImage>::InputImageType::Pointer
RieszImageFilter<TInputImage>
::ComputeRieszWeightedNormByEigenValues(const RieszComponentsImageType* rieszComponents, const typename EigenValuesImageType::Pointer eigenValues) const
{
  typename InputImageType::Pointer normOutput = InputImageType::New();
  normOutput->SetRegions(rieszComponents->GetBufferedRegion());
  normOutput->Allocate();
  normOutput->FillBuffer(NumericTraits<InputImagePixelType>::Zero);

  ImageRegionIterator<InputImageType> normOutputIt(
    normOutput, normOutput->GetRequestedRegion());
  ImageRegionConstIterator<RieszComponentsImageType> rieszIt(
    rieszComponents, rieszComponents->GetRequestedRegion());
  ImageRegionConstIterator<EigenValuesImageType> eigenValuesIt(
    eigenValues, eigenValues->GetRequestedRegion());

  InputImagePixelType localValue;
  itk::VariableLengthVector<InputImagePixelType> localEigenValues;
  localEigenValues.SetSize(ImageDimension);
  DirectionType weights;
  typename DirectionType::ValueType eigen_values_sum ;
  for(rieszIt.GoToBegin(), eigenValuesIt.GoToBegin(),
    normOutputIt.GoToBegin() ;
    !rieszIt.IsAtEnd() ;
    ++rieszIt, ++eigenValuesIt, ++normOutputIt )
    {
      eigen_values_sum = 0;
      localValue = 0;
      // Get the sum of eigenValues first
      for (unsigned int i = 0 ; i < ImageDimension ; ++i)
      {
        localEigenValues[i] = eigenValuesIt.Get()[i];
        eigen_values_sum += localEigenValues[i];
      }
      // std::cout << "Sum eigenValues: " << eigen_values_sum << std::endl;
      // Compute the weighted norm using eigenValues
      for (unsigned int i = 0 ; i < ImageDimension ; ++i)
      {
          weights[i] = localEigenValues[i] / eigen_values_sum;
          localValue += localValue + rieszIt.Get()[i] * rieszIt.Get()[i] *
            weights[i] * weights[i];
          // std::cout << "Weight[" << i << "]: " << weights[i] <<
          //   "localValue (temp)" << localValue << std::endl;
      }

      normOutputIt.Set( sqrt(localValue) );
    }

  return normOutput;
}
template< typename TInputImage>
std::pair<
  typename RieszImageFilter<TInputImage>::EigenVectorsImageType::Pointer,
  typename RieszImageFilter<TInputImage>::EigenValuesImageType::Pointer >
RieszImageFilter<TInputImage>
::ComputeEigenAnalysisMaximizingRieszComponents(const unsigned int & gaussian_window_radius, const float & gaussian_window_sigma, const RieszComponentsImageType* rieszComponents) const
{
  // For each pixel of eigenOut (size of TInput)
  //   For each RieszComponent:
  // Use NeighborhoodIterator in RieszComponents using gaussian_radius
  //  Weight value of pixels in the neighborhood using the GaussianImage
  //
  //  Compute Matrix (2D,Dimension x Dimension) eigenMatrix, which is the weighted contribution of the RieszComponents. J(x_0)[m][n] = Sum_each_neighbor_pixel_x(gaussian(x) * RieszComponent(x)[m] * RieszComponent(x)[n] )
  //   Compute eigenVectors and eigenValues of the matrix.
  //   Store them in the output.

  // Allocate outputs:

  typename EigenVectorsImageType::Pointer eigenVectorsImage = EigenVectorsImageType::New();
  typename EigenValuesImageType::Pointer eigenValuesImage = EigenValuesImageType::New();
  eigenVectorsImage->SetRegions(this->GetInput()->GetBufferedRegion());
  eigenVectorsImage->Allocate();
  eigenValuesImage->SetRegions(this->GetInput()->GetBufferedRegion());
  eigenValuesImage->Allocate();

  // Set GaussianImageSource
  Size<ImageDimension> radius;
  Size<ImageDimension> domainKernelSize;
  radius.Fill(gaussian_window_radius);
  domainKernelSize.Fill(2 * gaussian_window_radius + 1);
  typedef GaussianImageSource< InputImageType > GaussianSourceType;
  typename GaussianSourceType::Pointer gaussianImage =
    GaussianSourceType::New();
  typename GaussianSourceType::ArrayType mean;
  typename GaussianSourceType::ArrayType sigma;

  const typename RieszComponentsImageType::SpacingType inputSpacing =
    rieszComponents->GetSpacing();
  const typename RieszComponentsImageType::PointType inputOrigin =
    rieszComponents->GetOrigin();
  gaussianImage->SetSize( domainKernelSize );
  gaussianImage->SetSpacing(inputSpacing);
  gaussianImage->SetOrigin(inputOrigin);
  gaussianImage->SetScale(1.0);
  gaussianImage->SetNormalized(true);

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    mean[i] = inputSpacing[i] * radius[i] + inputOrigin[i]; // center pixel pos
    sigma[i] = gaussian_window_sigma;
    }
  gaussianImage->SetSigma(sigma);
  gaussianImage->SetMean(mean);
  gaussianImage->Update();
  std::cout << "GaussianImage" << std::endl;
  visualize::VisualizeITKImage(gaussianImage->GetOutput());
  // END Set GaussianImageSource

  ImageRegionIterator<EigenVectorsImageType> eigenVectorsIt(
      eigenVectorsImage,
      eigenVectorsImage->GetRequestedRegion());
  ImageRegionIterator<EigenValuesImageType> eigenValuesIt(
      eigenValuesImage,
      eigenValuesImage->GetRequestedRegion());
  ConstNeighborhoodIterator<RieszComponentsImageType> rieszIt(
      radius, rieszComponents,
      rieszComponents->GetRequestedRegion());
  ConstNeighborhoodIterator<typename GaussianSourceType::OutputImageType>
    gaussianIt(
      radius, gaussianImage->GetOutput(),
      gaussianImage->GetOutput()->GetRequestedRegion());

  typedef itk::Matrix<RealType, ImageDimension,ImageDimension> EigenMatrixType;
  typedef itk::FixedArray<RealType, ImageDimension> EigenValuesType;
  typedef itk::SymmetricEigenAnalysis<EigenMatrixType, EigenValuesType> SymmetricEigenAnalysisType;
  EigenMatrixType eigenMatrix;
  EigenMatrixType eigenVectors;
  EigenValuesType eigenValues;
  SymmetricEigenAnalysisType eigenSystem(ImageDimension);

  for(rieszIt.GoToBegin(), gaussianIt.GoToBegin(),
    eigenVectorsIt.GoToBegin(), eigenValuesIt.GoToBegin() ;
    !rieszIt.IsAtEnd() ;
    ++rieszIt, ++eigenVectorsIt, ++eigenValuesIt )
    {
      //Set the matrix
      eigenMatrix.Fill(0);
      for (unsigned int r = 0 ; r <= gaussian_window_radius ; ++r)
        for (unsigned int m = 0; m < ImageDimension ; ++m)
          for (unsigned int n = m ; n < ImageDimension ; ++n)
            if (r == 0)
              eigenMatrix[m][n] = eigenMatrix[n][m] +=
                gaussianIt.GetCenterPixel() +
                rieszIt.GetCenterPixel()[m] + rieszIt.GetCenterPixel()[n];
            else
              for (unsigned int axis = 0 ; axis < ImageDimension ; ++axis)
                {
                eigenMatrix[m][n] = eigenMatrix[n][m] +=
                  gaussianIt.GetNext(axis, r) +
                  rieszIt.GetNext(axis,r)[m] + rieszIt.GetNext(axis,r)[n];
                eigenMatrix[m][n] +=
                  gaussianIt.GetPrevious(axis, r) +
                  rieszIt.GetPrevious(r)[m] + rieszIt.GetPrevious(r)[n];
                }
      eigenSystem.ComputeEigenValuesAndVectors(eigenMatrix, eigenValues, eigenVectors );
      // Copy to Output
      eigenVectorsIt.Set(eigenVectors);
      eigenValuesIt.Set(eigenValues);

    } // end NeighborhoodIterator rieszIt

  return std::make_pair(eigenVectorsImage, eigenValuesImage);
  // return eigenOutputOriginal;
}

template< typename TInputImage>
typename RieszImageFilter<TInputImage>::RieszComponentsImageType::Pointer
RieszImageFilter<TInputImage>
::ComputeRieszComponentsWithMaximumResponse(const typename EigenVectorsImageType::Pointer eigenVectors, const RieszComponentsImageType* rieszComponents ) const
{
  typename RieszComponentsImageType::Pointer maxLocalRiesz =
    RieszComponentsImageType::New();
  maxLocalRiesz->SetRegions(rieszComponents->GetLargestPossibleRegion());
  maxLocalRiesz->SetVectorLength(ImageDimension);
  maxLocalRiesz->Allocate();

  ImageRegionIterator<RieszComponentsImageType> maxLocalRieszIt(
    maxLocalRiesz, maxLocalRiesz->GetRequestedRegion());
  ImageRegionConstIterator<RieszComponentsImageType> rieszIt(
    rieszComponents, rieszComponents->GetRequestedRegion());
  ImageRegionConstIterator<EigenVectorsImageType> eigenVectorsIt(
    eigenVectors, eigenVectors->GetRequestedRegion());

  itk::VariableLengthVector<InputImagePixelType> rieszProjection;
  rieszProjection.SetSize(ImageDimension);
  DirectionType direction;
  for(rieszIt.GoToBegin(), eigenVectorsIt.GoToBegin(),
    maxLocalRieszIt.GoToBegin() ;
    !rieszIt.IsAtEnd() ;
    ++rieszIt, ++eigenVectorsIt, ++maxLocalRieszIt )
    {
      for (unsigned int column = 0 ; column < ImageDimension ; ++column){
        for (unsigned int i = 0 ; i < ImageDimension ; ++i)
          direction[i] = eigenVectorsIt.Get()[i][column];
        rieszProjection[column] = ComputeLocalRieszProjection(direction, rieszIt);
      }
      maxLocalRieszIt.Set( rieszProjection );
    }

  return maxLocalRiesz;
}

template< typename TInputImage>
typename TInputImage::Pointer
RieszImageFilter<TInputImage>
::ComputeLocalPhaseInDirection(const DirectionType & unitary_direction,
    const InputImageType* rieszReal, const RieszComponentsImageType* rieszComponents) const
{
  typename InputImageType::Pointer rieszProjection = ComputeRieszProjection(unitary_direction, rieszComponents);
  rieszProjection->DisconnectPipeline();
  typedef itk::DivideImageFilter<TInputImage,TInputImage, TInputImage> DivideFilterType;
  typename DivideFilterType::Pointer divideFilter = DivideFilterType::New();
  divideFilter->SetInput1(rieszProjection);
  divideFilter->SetInput2(rieszReal);
  divideFilter->Update();

  typedef itk::AtanImageFilter<TInputImage,TInputImage> AtanFilterType;
  typename AtanFilterType::Pointer atanFilter = AtanFilterType::New();
  atanFilter->SetInput(divideFilter->GetOutput());
  atanFilter->Update();
  typename InputImageType::Pointer output = atanFilter->GetOutput();
  output->DisconnectPipeline();
  return output;
}

template< typename TInputImage>
typename TInputImage::Pointer
RieszImageFilter<TInputImage>
::ComputeRieszProjection(const DirectionType & unitary_direction,
    const RieszComponentsImageType* rieszComponents) const
{
  typedef itk::ImageDuplicator< InputImageType > DuplicatorType;
  typename DuplicatorType::Pointer duplicator= DuplicatorType::New();
  duplicator->SetInputImage(this->GetInput());
  duplicator->Update();
  typename InputImageType::Pointer output = duplicator->GetModifiableOutput();

  typedef VectorIndexSelectionCastImageFilter<RieszComponentsImageType,InputImageType> CastIndexType;
  typename CastIndexType::Pointer castIndex = CastIndexType::New();
  castIndex->SetInput(rieszComponents);

  typedef itk::MultiplyImageFilter<TInputImage,TInputImage,TInputImage> MultiplyImageFilterType;
  typename MultiplyImageFilterType::Pointer multiplyByConstant = MultiplyImageFilterType::New();

  typedef AddImageFilter<InputImageType,InputImageType,InputImageType> AddImageFilterType;
  typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();

  for (unsigned int i = 0 ; i<ImageDimension ; ++i){
    castIndex->SetIndex(i);
    castIndex->Update();
    multiplyByConstant->SetInput(castIndex->GetOutput());
    multiplyByConstant->SetConstant(unitary_direction[i]);
    multiplyByConstant->Update();

    addFilter->SetInput1(output);
    addFilter->SetInput2(multiplyByConstant->GetOutput());
    addFilter->Update();
    output = addFilter->GetOutput();
    output->DisconnectPipeline();
  }
  return output;

}
template< typename TInputImage>
typename RieszImageFilter<TInputImage>::RealType
RieszImageFilter<TInputImage>
::ComputeLocalRieszProjection(const DirectionType & direction, const ImageConstIterator<RieszComponentsImageType> & rieszIt ) const
{
  RealType rieszProjection = 0;
  for (unsigned int i = 0 ; i < ImageDimension ; ++i)
    rieszProjection += direction[i] * rieszIt.Get()[i] ;
  return rieszProjection;
}

template< typename TInputImage >
void
RieszImageFilter< TInputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Sigma of Gaussian Derivative : " << m_SigmaGaussianDerivative
     << std::endl;
  os << indent << "Size Square: " << m_inputSizeSquare <<  std::endl;
  os << indent << "Spacing Square: " << m_inputSpacingSquare << std::endl;
}

template< typename TInputImage>
DataObject::Pointer RieszImageFilter<TInputImage>
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
    default:
      std::cerr << "No output " << idx << std::endl;
      output = NULL;
      break;
    }
  return output.GetPointer();
}

template< class TInputImage>
TInputImage* RieszImageFilter<TInputImage>::GetOutputReal()
{
  return dynamic_cast< TInputImage * >(
           this->ProcessObject::GetOutput(0) );
}

template< class TInputImage>
typename RieszImageFilter<TInputImage>::RieszComponentsImageType* RieszImageFilter<TInputImage>::GetOutputRieszComponents()
{
  return dynamic_cast< RieszComponentsImageType * >(
           this->ProcessObject::GetOutput(1) );
}

template< class TInputImage>
TInputImage* RieszImageFilter<TInputImage>::GetOutputRieszNorm()
{
  return dynamic_cast< TInputImage * >(
           this->ProcessObject::GetOutput(2) );
}

} // end namespace itk
#endif
