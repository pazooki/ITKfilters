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
#ifndef itkFrequencyExpandImageFilter_hxx
#define itkFrequencyExpandImageFilter_hxx

#include "itkFrequencyExpandImageFilter.h"
#include "itkObjectFactory.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"
#include "Ind2Sub.h"
#include "itkPasteImageFilter.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{
/**
 * Default constructor
 */
template< typename TInputImage, typename TOutputImage >
FrequencyExpandImageFilter< TInputImage, TOutputImage >
::FrequencyExpandImageFilter()
{
  // Set default factors to 1
  for ( unsigned int j = 0; j < ImageDimension; j++ )
    {
    m_ExpandFactors[j] = 1;
    }
}

/**
 * Standard "PrintSelf" method
 */
template< typename TInputImage, typename TOutputImage >
void
FrequencyExpandImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  unsigned int j;
  os << indent << "ExpandFactors: [";
  for ( j = 0; j < ImageDimension - 1; j++ )
    {
    os << m_ExpandFactors[j] << ", ";
    }
  os << m_ExpandFactors[j] << "]" << std::endl;
}

/**
 * Set expand factors from a single unsigned int
 */
template< typename TInputImage, typename TOutputImage >
void
FrequencyExpandImageFilter< TInputImage, TOutputImage >
::SetExpandFactors(
  const unsigned int factor)
{
  unsigned int j;

  for ( j = 0; j < ImageDimension; j++ )
    {
    if ( factor != m_ExpandFactors[j] ) { break; }
    }
  if ( j < ImageDimension )
    {
    this->Modified();
    for ( j = 0; j < ImageDimension; j++ )
      {
      m_ExpandFactors[j] = factor;
      if ( m_ExpandFactors[j] < 1 ) { m_ExpandFactors[j] = 1; }
      }
    }
}

template< typename TInputImage, typename TOutputImage >
void
FrequencyExpandImageFilter< TInputImage, TOutputImage >
::GenerateData()
{
  const InputImageType * inputPtr  = this->GetInput();
  OutputImagePointer     outputPtr = this->GetOutput();

  typename TInputImage::SizeType inputSize = inputPtr->GetLargestPossibleRegion().GetSize();
  typename TOutputImage::SizeType outputSize = outputPtr->GetLargestPossibleRegion().GetSize();

  typename TInputImage::SizeType halfInputSize;
  for ( unsigned int i=0; i < TInputImage::ImageDimension; ++i )
  {
    // TODO Do we want to check shrinkfactor is integer divisor of inputSize?
    // YES, only work for even inputImg
    halfInputSize[i] = inputSize[i] / 2 ;
  }
  const typename TOutputImage::IndexType indexRequested = outputPtr->GetLargestPossibleRegion().GetIndex();

  // Manage ImageDimension array linearly:{{{
  std::array<unsigned int , ImageDimension> nsizes;
  unsigned int numberOfRegions = 1;
  for (unsigned int dim = 0 ; dim < ImageDimension ; ++dim)
    {
    nsizes[dim] = 2;
    numberOfRegions*=nsizes[dim];
    }
  std::array<unsigned int, ImageDimension> subIndices;
  /// }}}


  // Prepare filter to paste the different regions into output.
  typedef itk::PasteImageFilter<OutputImageType> PasteFilterType;
  typename PasteFilterType::Pointer pasteFilter = PasteFilterType::New();
  pasteFilter->SetSourceImage(inputPtr);
  pasteFilter->SetDestinationImage(outputPtr);
  pasteFilter->InPlaceOn();

  typedef typename OutputImageType::RegionType RegionType;
  ProgressReporter progress(this, 0, numberOfRegions );

  for (unsigned int n = 0 ; n < numberOfRegions ; ++n)
    {
    subIndices = Ind2Sub<ImageDimension>(n, nsizes);
    RegionType zoneRegion;
    typename OutputImageType::SizeType zoneSize;
    typename TInputImage::IndexType  inputIndex  = indexRequested;
    typename TOutputImage::IndexType outputIndex = indexRequested;
    for (unsigned int dim = 0 ; dim < ImageDimension ; ++dim)
      {
      if(subIndices[dim] == 1) // TODO, only for even now.
        {
        // Negative frequencies
        inputIndex[dim] = halfInputSize[dim] ;
        zoneSize[dim] = halfInputSize[dim] ;
        outputIndex[dim] = outputSize[dim] - halfInputSize[dim] ;
        }
      else
        {
        inputIndex[dim] = 0 ;
        zoneSize[dim] = halfInputSize[dim] ;
        outputIndex[dim] = 0 ;
        }
      }
    zoneRegion.SetIndex(inputIndex);
    zoneRegion.SetSize(zoneSize);
    itkDebugMacro( << "n:" << n << " region: " << zoneRegion);

    pasteFilter->SetSourceRegion(zoneRegion);
    pasteFilter->SetDestinationIndex(outputIndex);
    if (n == numberOfRegions - 1) // Graft the output.
      {
      pasteFilter->GraftOutput(outputPtr);
      pasteFilter->Update();
      this->GraftOutput(pasteFilter->GetOutput());
      }
    else // update output
      {
      pasteFilter->Update();
      outputPtr = pasteFilter->GetOutput();
      }
    progress.CompletedPixel();
    }
  /* TODO: Do you have to "recalculate" DC component zero frequency bins?
  *           modifyFFT --> What is the new DC?
  * ANSWER: I don't think you can calculate modified DC from FFT.DC components
  * can be set to zero BEFORE performing FFT.
  * data = data - mean(data) ---> FFT ---> DC = 0
  */
  // Set zero value to DC: mean real value
  // itk::ImageRegionConstIterator<OutputImageType> outIt(
  //   this->GetOutput(),
  //   this->GetOutput()->GetLargestPossibleRegion());
  // outIt.GoToBegin();
  // size_t linearOutputSize(1);
  // for ( unsigned int i=0; i < TInputImage::ImageDimension; ++i )
  //   linearOutputSize *= outputSize[i];
  // typename OutputImageType::PixelType mean(0);
  // // Skip from the mean calculation the first index that stores the old zero frequency value.
  // ++outIt;
  // while(!outIt.IsAtEnd())
  //   {
  //   // mean += outIt.Get() * (1.0/linearOutputSize);
  //   mean += outIt.Get();
  //   ++outIt;
  //   }
  // this->GetOutput()->SetPixel(this->GetOutput()->GetLargestPossibleRegion().GetIndex(), mean);

}

/**
 * GenerateInputRequesteRegion
 */
template< typename TInputImage, typename TOutputImage >
void
FrequencyExpandImageFilter< TInputImage, TOutputImage >
::GenerateInputRequestedRegion()
{
  // Call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // Get pointers to the input and output
  InputImagePointer inputPtr =
    const_cast< TInputImage * >( this->GetInput() );
  OutputImagePointer outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // We need to compute the input requested region (size and start index)
  unsigned int i;
  const typename TOutputImage::SizeType & outputRequestedRegionSize =
    outputPtr->GetRequestedRegion().GetSize();
  const typename TOutputImage::IndexType & outputRequestedRegionStartIndex =
    outputPtr->GetRequestedRegion().GetIndex();

  typename TInputImage::SizeType inputRequestedRegionSize;
  typename TInputImage::IndexType inputRequestedRegionStartIndex;

  /**
   * inputRequestedSize = (outputRequestedSize / ExpandFactor) + 1)
   * The extra 1 above is to take care of edge effects when streaming.
   */
  for ( i = 0; i < TInputImage::ImageDimension; i++ )
    {
    inputRequestedRegionSize[i] =
      (SizeValueType)std::ceil( (double)outputRequestedRegionSize[i]
                      / (double)m_ExpandFactors[i] ) + 1;

    inputRequestedRegionStartIndex[i] =
      (SizeValueType)std::floor( (double)outputRequestedRegionStartIndex[i]
                       / (double)m_ExpandFactors[i] );
    }

  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion.SetSize(inputRequestedRegionSize);
  inputRequestedRegion.SetIndex(inputRequestedRegionStartIndex);

  // Make sure the requested region is within largest possible.
  inputRequestedRegion.Crop( inputPtr->GetLargestPossibleRegion() );

  // Set the input requested region.
  inputPtr->SetRequestedRegion(inputRequestedRegion);
}

/**
 * GenerateOutputInformation
 */
template< typename TInputImage, typename TOutputImage >
void
FrequencyExpandImageFilter< TInputImage, TOutputImage >
::GenerateOutputInformation()
{
  // Call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();

  // Get pointers to the input and output
  InputImagePointer inputPtr =
    const_cast< TInputImage * >( this->GetInput() );
  OutputImagePointer outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // We need to compute the output spacing, the output image size, and the
  // output image start index
  const typename TInputImage::SpacingType &
  inputSpacing = inputPtr->GetSpacing();
  const typename TInputImage::SizeType &   inputSize =
    inputPtr->GetLargestPossibleRegion().GetSize();
  const typename TInputImage::IndexType &  inputStartIndex =
    inputPtr->GetLargestPossibleRegion().GetIndex();
  const typename TInputImage::PointType &
  inputOrigin = inputPtr->GetOrigin();

  typename TOutputImage::SpacingType outputSpacing;
  typename TOutputImage::SizeType outputSize;
  typename TOutputImage::IndexType outputStartIndex;
  typename TOutputImage::PointType outputOrigin;

  typename TInputImage::SpacingType inputOriginShift;

  for ( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
    {
    outputSpacing[i] = inputSpacing[i] / (float)m_ExpandFactors[i];
    outputSize[i] = inputSize[i] * (SizeValueType)m_ExpandFactors[i];
    outputStartIndex[i] = inputStartIndex[i] * (IndexValueType)m_ExpandFactors[i];
    const double fraction = (double)( m_ExpandFactors[i] - 1 ) / (double)m_ExpandFactors[i];
    inputOriginShift[i] = -( inputSpacing[i] / 2.0 ) * fraction;
    }

  const typename TInputImage::DirectionType inputDirection = inputPtr->GetDirection();
  const typename TOutputImage::SpacingType outputOriginShift = inputDirection * inputOriginShift;

  outputOrigin = inputOrigin + outputOriginShift;

  outputPtr->SetSpacing(outputSpacing);
  outputPtr->SetOrigin(outputOrigin);

  typename TOutputImage::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize(outputSize);
  outputLargestPossibleRegion.SetIndex(outputStartIndex);

  outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);
}
} // end namespace itk

#endif
