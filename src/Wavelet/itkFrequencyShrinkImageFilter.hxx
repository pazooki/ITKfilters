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
#ifndef itkFrequencyShrinkImageFilter_hxx
#define itkFrequencyShrinkImageFilter_hxx

#include "itkFrequencyShrinkImageFilter.h"
#include "itkImageScanlineIterator.h"
#include "itkProgressReporter.h"
#include <numeric>
#include <functional>
#include <Ind2Sub.h>
#include <itkPasteImageFilter.h>

namespace itk
{

template< class TInputImage, class TOutputImage >
FrequencyShrinkImageFilter< TInputImage, TOutputImage >
::FrequencyShrinkImageFilter()
{
  for ( unsigned int j = 0; j < ImageDimension; j++ )
    {
    m_ShrinkFactors[j] = 1;
    }
}

template< class TInputImage, class TOutputImage >
void
FrequencyShrinkImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Shrink Factor: ";
  for ( unsigned int j = 0; j < ImageDimension; j++ )
    {
    os << m_ShrinkFactors[j] << " ";
    }
  os << std::endl;
}

template< class TInputImage, class TOutputImage >
void
FrequencyShrinkImageFilter< TInputImage, TOutputImage >
::SetShrinkFactors(unsigned int factor)
{
  unsigned int j;

  for ( j = 0; j < ImageDimension; j++ )
    {
    if ( factor != m_ShrinkFactors[j] ) { break; }
    }
  if ( j < ImageDimension )
    {
    this->Modified();
    for ( j = 0; j < ImageDimension; j++ )
      {
      m_ShrinkFactors[j] = factor;
      if ( m_ShrinkFactors[j] < 1 )
        {
        m_ShrinkFactors[j] = 1;
        }
      }
    }
}

template< class TInputImage, class TOutputImage >
void
FrequencyShrinkImageFilter< TInputImage, TOutputImage >
::SetShrinkFactor(unsigned int i, unsigned int factor)
{
  if ( m_ShrinkFactors[i] == factor )
    {
    return;
    }

  this->Modified();
  m_ShrinkFactors[i] = factor;
}

template <class TInputImage, class TOutputImage>
void
FrequencyShrinkImageFilter<TInputImage,TOutputImage>
::GenerateData()
{

  this->AllocateOutputs();
  // Get the input and output pointers
  const InputImageType * inputPtr = this->GetInput();
  typename OutputImageType::Pointer      outputPtr = this->GetOutput();
  outputPtr->SetBufferedRegion(outputPtr->GetLargestPossibleRegion());

  typename TInputImage::SizeType inputSize = inputPtr->GetLargestPossibleRegion().GetSize();
  typename TOutputImage::SizeType outputSize = outputPtr->GetLargestPossibleRegion().GetSize();
  typename TInputImage::SizeType lowPartSize;
  for ( unsigned int i=0; i < TInputImage::ImageDimension; ++i )
  {
    // TODO Do we want to check shrinkfactor is integer divisor of inputSize?
    // YES, only work for even inputImg
    lowPartSize[i] = outputSize[i] / 2 ;
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
        inputIndex[dim] = inputSize[dim] - lowPartSize[dim] ;
        zoneSize[dim] = lowPartSize[dim] ;
        outputIndex[dim] = lowPartSize[dim] ;
        }
      else
        {
        inputIndex[dim] = 0 ;
        zoneSize[dim] = lowPartSize[dim] ;
        outputIndex[dim] = 0;
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
}

template <class TInputImage, class TOutputImage>
void
FrequencyShrinkImageFilter<TInputImage,TOutputImage>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  InputImagePointer  inputPtr = const_cast<TInputImage *> (this->GetInput() );
  OutputImagePointer outputPtr = this->GetOutput();

  // The filter chops high frequencys [0 1...H,H-1 H-2...1], and averagage low frequencies. We need the whole input image, indepently of the RequestedRegion.
  inputPtr->SetRequestedRegion( inputPtr->GetLargestPossibleRegion() );
}

template <class TInputImage, class TOutputImage>
void
FrequencyShrinkImageFilter<TInputImage,TOutputImage>
::GenerateOutputInformation()
{
  // Call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();

  // Get pointers to the input and output
  InputImageConstPointer inputPtr  = this->GetInput();
  OutputImagePointer     outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // Compute the output spacing, the output image size, and the
  // output image start index
  const typename TInputImage::SpacingType & inputSpacing = inputPtr->GetSpacing();
  const typename TInputImage::SizeType &   inputSize = inputPtr->GetLargestPossibleRegion().GetSize();
  const typename TInputImage::IndexType &  inputStartIndex = inputPtr->GetLargestPossibleRegion().GetIndex();

  ContinuousIndex<double,ImageDimension> inputIndexOutputOrigin;

  typename TOutputImage::SpacingType outputSpacing(inputSpacing);
  typename TOutputImage::SizeType outputSize;
  typename TOutputImage::PointType outputOrigin;
  typename TOutputImage::IndexType outputStartIndex;

  for ( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
    {
    outputSpacing[i] *= m_ShrinkFactors[i];

    inputIndexOutputOrigin[i] = 0.5*(m_ShrinkFactors[i]-1);

    outputStartIndex[i] = Math::Ceil<SizeValueType>(inputStartIndex[i]/static_cast<double>( m_ShrinkFactors[i]) );

    // Round down so that all output pixels fit input input region
    outputSize[i] =
      Math::Floor<SizeValueType>( (double)(inputSize[i] - outputStartIndex[i]*m_ShrinkFactors[i]+inputStartIndex[i])
                                  / (double)m_ShrinkFactors[i]);

    if ( outputSize[i] < 1 )
      {
      itkExceptionMacro("InputImage is too small! An output pixel does not map to a whole input bin.");
      }

    }

  inputPtr->TransformContinuousIndexToPhysicalPoint(inputIndexOutputOrigin, outputOrigin);

  outputPtr->SetSpacing(outputSpacing);
  outputPtr->SetOrigin(outputOrigin);

  // Set region
  typename TOutputImage::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize(outputSize);
  outputLargestPossibleRegion.SetIndex(outputStartIndex);

  outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);

}

} // end namespace itk

#endif
