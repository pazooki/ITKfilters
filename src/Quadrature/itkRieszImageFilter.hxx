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

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRieszImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkForwardFFTImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkDiscreteGaussianDerivativeImageFunction.h"
#include "itkMultiplyImageFilter.h"
#include <complex>


namespace itk
{
template< typename TInputImage, typename TOutputImage >
RieszImageFilter< TInputImage, TOutputImage >
::RieszImageFilter()
{
  // m_Height = 2;
  // m_NumberOfIterationsUsed = 1;
  // m_FullyConnected = false;
}
template< typename TInputImage, typename TOutputImage >
void
RieszImageFilter< TInputImage, TOutputImage >
::GenerateData()
{
  // this->AllocateOutputs();
  // Allocate the output
  OutputImageType* outputPtr = this->GetOutput();
  // allocate the output buffer
  outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
  outputPtr->Allocate();

  std::cout << "FFT forward" << std::endl;
  typedef ForwardFFTImageFilter < InputImageType> FFTFilterType;
  typename FFTFilterType::Pointer fftFilter = FFTFilterType::New();
  fftFilter->SetInput( this->GetInput() );
  fftFilter->Update();

  // variance of gaussian derivative-> TODO put in member.
  double m_variance_gaussian_derivative = 5.0;

  ImageRegionIteratorWithIndex< OutputImageType >
  outIt( outputPtr, outputPtr->GetRequestedRegion() );
  // Walk the output requested region
  double w_mod = 0;
  for ( outIt.GoToBegin() ; !outIt.IsAtEnd(); ++outIt )
    {
    const typename OutputImageType::IndexType index = outIt.GetIndex();
    // The position at which the function is evaluated
    typename OutputImageType::PointType evalPoint; // w = (u,v,w);
    outputPtr->TransformIndexToPhysicalPoint(index, evalPoint);
    w_mod = sqrt( evalPoint[0]*evalPoint[0] +
                  evalPoint[1]*evalPoint[1] +
                  evalPoint[2]*evalPoint[2] );

    // Complex value:
    // const std::complex<InputImagePixelType> value =
    const OutputImagePixelType value = outIt.Get() *
      static_cast<InputImagePixelType>(w_mod * std::exp(m_variance_gaussian_derivative * w_mod) );

    // Set the pixel value to the function value
    outIt.Set( value );
    }


  // ---- start bandpass filter ----
  // Gaussian Isotropic Filter in frequency (or Gabor is another option)
    // typedef DiscreteGaussianDerivativeImageFunction<OutputImageType, OutputImageType> DiscreteGaussianDerivativeType;
    // typename DiscreteGaussianDerivativeType::Pointer discreteGaussian = DiscreteGaussianDerivativeType::New();
    // discreteGaussian->SetVariance(m_variance);
    // discreteGaussian->SetUseImageSpacing(false); // Variance in pixels.
    // discreteGaussian->SetOrder(1); // Derivative order.
    // discreteGaussian->SetNormalizeAcrossScale( false );
    // discreteGaussian->SetInput(fftFilter->GetOutput());
    // discreteGaussian->Update();

    // typedef MultiplyImageFilter <OutputImageType,InputImageType> MultiplyFilterType;
    // typename MultiplyFilterType::Pointer multiplier =
    //   MultiplyFilterType::New();
    // multiplier->SetInput1( fftFilter->GetOutput() );
    // multiplier->SetInput2( discreteGaussian->GetOutput() );
    // multiplier->Update();
    // ---- end bandpass filter ---/

  // // construct a marker image to manipulate using reconstruction by
  // // dilation. the marker image is the input image minus the height
  // // parameter.
  // typedef ShiftScaleImageFilter< TInputImage, TInputImage > ShiftFilterType;
  // typename ShiftFilterType::Pointer shift = ShiftFilterType::New();
  // shift->SetInput( this->GetInput() );
  // shift->SetShift( -1.0 * static_cast< typename ShiftFilterType::RealType >( m_Height ) );
  //
  // // Delegate to a geodesic dilation filter.
  // //
  // //
  // typename ReconstructionByDilationImageFilter< TInputImage, TInputImage >::Pointer
  // dilate =
  //   ReconstructionByDilationImageFilter< TInputImage, TInputImage >::New();
  //
  // // Create a process accumulator for tracking the progress of this minipipeline
  // ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  // progress->SetMiniPipelineFilter(this);
  // progress->RegisterInternalFilter(dilate, 1.0f);
  //
  // // set up the dilate filter
  // //dilate->RunOneIterationOff();             // run to convergence
  // dilate->SetMarkerImage( shift->GetOutput() );
  // dilate->SetMaskImage( this->GetInput() );
  // dilate->SetFullyConnected(m_FullyConnected);
  //
  // Must cast to the output type
  // typename CastImageFilter< TInputImage, TOutputImage >::Pointer cast =
  //   CastImageFilter< TInputImage, TOutputImage >::New();
  // cast->SetInput( dilate->GetOutput() );
  // cast->InPlaceOn();

  // graft our output to the cast filter to force the proper regions
  // to be generated
  // cast->GraftOutput( this->GetOutput() );

  // reconstruction by dilation
  // cast->Update();

  // graft the output of the dilate filter back onto this filter's
  // output. this is needed to get the appropriate regions passed
  // back.
  // this->GraftOutput( cast->GetOutput() );
}

template< typename TInputImage, typename TOutputImage >
void
RieszImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  // os << indent << "Height of local maxima (contrast): "
  //    << static_cast< typename NumericTraits< InputImagePixelType >::PrintType >( m_Height )
  //    << std::endl;
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
      output = ( InputImageType::New() ).GetPointer();
      break;
    case 2:
      output = ( InputImageType::New() ).GetPointer();
      break;
    case 3:
      output = ( InputImageType::New() ).GetPointer();
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
