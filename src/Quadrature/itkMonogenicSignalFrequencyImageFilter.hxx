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
#ifndef itkMonogenicSignalFrequencyImageFilter_hxx
#define itkMonogenicSignalFrequencyImageFilter_hxx
#include "itkMonogenicSignalFrequencyImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkRieszFrequencyFunction.h"
namespace itk
{
template< typename TInputImage
          typename TFrequencyRegionIterator>
MonogenicSignalFrequencyImageFilter<TInputImage, TFrequencyRegionIterator>
::MonogenicSignalFrequencyImageFilter()
{
}

template< typename TInputImage
          typename TFrequencyRegionIterator>
MonogenicSignalFrequencyImageFilter<TInputImage, TFrequencyRegionIterator>
::GenerateOutputInformation(void)
{
  this->Superclass::GenerateOutputInformation();

  OutputImageType *output = this->GetOutput();
  output->SetNumberOfComponentsPerPixel( ImageDimension + 1 );
}

template< typename TInputImage >
void
MonogenicSignalFrequencyImageFilter< TInputImage >
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                       ThreadIdType threadId)
{
  // Allocate the outputs
  this->AllocateOutputs();

  InputFrequencyImageRegionIterator   inFreqIt(this->GetInput() , outputRegionForThread);
  ImageRegionIterator< OutputImageType > outIt(this->GetOutput(), outputRegionForThread);

  typedef itk::RieszFrequencyFunction<FloatType, ImageDimension> RieszFunctionType;
  typename RieszFunctionType::Pointer evaluator = RieszFunctionType::New();

  for (inFreqIt.GoToBegin(), outIt.GoToBegin(); !inFreqIt.IsAtEnd(); ++inFreqIt, ++outIt)
    {
    typename RieszFunctionType::OutputArrayType evaluatedArray =
      evaluator->EvaluateArray(inFreqIt.GetFrequencyIndex());
    typename OutputImageType::PixelType out_value = outIt.Get();
    out_value[0] = inFreqIt.Get();
    for (unsigned int dir = 0; dir < ImageDimension; ++dir)
      {
      // This is a complex number multiplication.
      value[dir + 1] = inFreqIt.Get() * evaluatedArray[dir];
      }
    outIt.Set(value);
    }
}
} // end namespace itk
#endif
