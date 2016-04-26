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
// #include "visualize_functions.h" // TODO REMOVE
#include <array>
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRieszImageFilter.h"
#include <itkComposeImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
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
#include <itkGaussianOperator.h>

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
  evenPart->SetBufferedRegion( evenPart->GetLargestPossibleRegion() );
  evenPart->Allocate();
  RieszComponentsImageType* rieszComponents = this->GetOutputRieszComponents();
  rieszComponents->SetBufferedRegion( rieszComponents->GetRequestedRegion() );
  rieszComponents->Allocate();
  InputImageType* normPart = this->GetOutputRieszNorm();
  normPart->SetBufferedRegion( normPart->GetRequestedRegion() );
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
::ComputeRieszNorm(const RieszComponentsImageType* rieszComponents) const
{
  // typedef itk::ImageDuplicator< InputImageType > DuplicatorType;
  // typename DuplicatorType::Pointer duplicator= DuplicatorType::New();
  // duplicator->SetInputImage(this->GetInput());
  // duplicator->Update();
  typename InputImageType::Pointer output = InputImageType::New();
  output->SetRegions(rieszComponents->GetBufferedRegion());
  output->Allocate();
  output->FillBuffer(NumericTraits<InputImagePixelType>::Zero);
  std::cout << "End Duplicate" << std::endl;

  typedef VectorIndexSelectionCastImageFilter<RieszComponentsImageType,InputImageType> CastIndexType;
  typename CastIndexType::Pointer castIndex = CastIndexType::New();
  castIndex->SetInput(rieszComponents);

  typedef SquareImageFilter<InputImageType,InputImageType> SquareImageFilterType;
  typename SquareImageFilterType::Pointer squareFilter = SquareImageFilterType::New();

  typedef AddImageFilter<InputImageType,InputImageType,InputImageType> AddImageFilterType;
  typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();

  for (unsigned int i = 0 ; i<ImageDimension ; ++i){
    castIndex->SetIndex(i);
    castIndex->Update();
    squareFilter->SetInput(castIndex->GetOutput());
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
typename RieszImageFilter<TInputImage>::EigenImageType::Pointer
RieszImageFilter<TInputImage>
::ComputeEigenVectorsMaximizingRieszComponents(const unsigned int & gaussian_window_radius, const RieszComponentsImageType* rieszComponents) const
{
  typename EigenImageType::Pointer eigenOutputOriginal = EigenImageType::New();
 // typename EigenImageType::Region region;
  // eigenOutputOriginal->SetBufferedRegion( rieszComponents->GetRequestedRegion() );
  eigenOutputOriginal->Allocate();


  std::array<InputImagePointer, ImageDimension> comps;
  typedef VectorIndexSelectionCastImageFilter<RieszComponentsImageType,InputImageType> CastIndexType;
  typename CastIndexType::Pointer castIndex = CastIndexType::New();
  castIndex->SetInput(rieszComponents);
  for (unsigned int i = 0 ; i < ImageDimension ; ++i)
  {
    castIndex->SetIndex(i);
    comps[i] = castIndex->GetOutput();
    comps[i]->DisconnectPipeline();
  }

  // typename RieszComponentsArrayImageType::Pointer rieszArray = RieszComponentsArrayImageType::New();
  // rieszArray->SetBufferedRegion(rieszComponents->GetBufferedRegion());
  // rieszArray->Allocate();
  // ImageRegionIterator<RieszComponentsArrayImageType> rieszArrayIt(
  //     rieszArray, rieszArray->GetRequestedRegion());
  // ImageRegionConstIterator<RieszComponentsImageType> rieszIt(
  //     rieszComponents, rieszComponents->GetRequestedRegion());
  // for(rieszArrayIt.GoToBegin(), rieszIt.GoToBegin() ;
  //     !rieszArrayIt.IsAtEnd();
  //     ++rieszArrayIt, ++rieszIt)
  // {
  //   auto g = rieszIt.Get();
  //   double carray[3] =  {g[0], g[1], g[2] } ;
  //   rieszArrayIt.Set(carray);
  // }
  // For each pixel of eigenOut (size of TInput)
  //   For each RieszComponent:
  //    Get neighborhood of each pixel of size equal to Gaussian Window.
  //    Convolve neighborhood with the window.
  //   Compute Matrix J (2D,Dimension x Dimension) using all convolved neighborhoods.
  //   Compute eigenVectors and eigenValues of matrix.
  //   Store them in eigenOutput.

  // TODO look ImageVectorToImageAdaptor for using the GaussianOperator.
  NeighborhoodInnerProduct<InputImageType> innerProduct;
  GaussianOperator<typename InputImageType::PixelType, ImageDimension> gaussianOperator;
  Size<ImageDimension> radius;
  radius.Fill(gaussian_window_radius);
  gaussianOperator.CreateToRadius(radius);

  ConstNeighborhoodIterator<RieszComponentsImageType> rieszNIt(
      gaussianOperator.GetRadius(), rieszComponents,
      rieszComponents->GetRequestedRegion());
  ImageRegionIterator<EigenImageType> eigenIt(
      eigenOutputOriginal, eigenOutputOriginal->GetRequestedRegion());

  // NeighborhoodInnerProduct<InputImageType> innerProduct;
  // GaussianOperator<typename InputImageType::PixelType, ImageDimension> gaussianOperator;
  // Size<ImageDimension> radius;
  // radius.Fill(gaussian_window_radius);
  // gaussianOperator.CreateToRadius(radius);
  //
  // ImageRegionIterator<EigenImageType> eigenIt(
  //     eigenOutputOriginal, eigenOutputOriginal->GetRequestedRegion());

  // ConstNeighborhoodIterator<RieszComponentsImageType> rieszIt(
  //     gaussianOperator.GetRadius(), rieszComponents,
  //     rieszComponents->GetRequestedRegion());

  while(!eigenIt.IsAtEnd())
  {

    ++eigenIt;
  }
  // for(rieszIt.GoToBegin(), eigenIt.GoToBegin() ;
  //     !rieszIt.IsAtEnd() ;
  //     ++rieszIt , ++eigenIt)
  // {
  //   auto a = innerProduct(rieszIt, gaussianOperator);
  //
  //   eigenIt.Set();
  // } // end NeighborhoodIterator rieszIt

  return eigenOutputOriginal;
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
::MakeOutput(unsigned int idx)
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
