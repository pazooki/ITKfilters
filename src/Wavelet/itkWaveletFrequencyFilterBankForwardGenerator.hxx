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
#ifndef itkWaveletFrequencyFilterBankForwardGenerator_hxx
#define itkWaveletFrequencyFilterBankForwardGenerator_hxx
#include "itkWaveletFrequencyFilterBankForwardGenerator.h"
#include "itkNumericTraits.h"

namespace itk {
template < typename TOutputImage, typename TWaveletFunction>
WaveletFrequencyFilterBankForwardGenerator< TOutputImage, TWaveletFunction>
::WaveletFrequencyFilterBankForwardGenerator()
  : m_HighPassSubBands(0),
    m_ScaleFactor(1)
{
  this->SetHighPassSubBands(1);
}


template < typename TOutputImage, typename TWaveletFunction>
void WaveletFrequencyFilterBankForwardGenerator< TOutputImage, TWaveletFunction>
::SetHighPassSubBands(unsigned int k)
{
  if ( m_HighPassSubBands == k ) return;

  this->m_HighPassSubBands = k;
  this->SetNumberOfRequiredOutputs(k + 1);
  this->Modified();

  for (unsigned int band = 0; band < this->m_HighPassSubBands + 1; ++band)
    this->SetNthOutput(band, this->MakeOutput(band));
}

template < typename TOutputImage, typename TWaveletFunction>
void WaveletFrequencyFilterBankForwardGenerator< TOutputImage, TWaveletFunction>
::PrintSelf(std::ostream &os, Indent indent) const {
  Superclass::PrintSelf(os, indent);

  os << indent << "HighPassSubBands: " << this->m_HighPassSubBands
     << std::endl;
}

/* ******* Get Outputs *****/
template < typename TOutputImage, typename TWaveletFunction>
typename WaveletFrequencyFilterBankForwardGenerator< TOutputImage, TWaveletFunction>
::OutputImagePointer
WaveletFrequencyFilterBankForwardGenerator< TOutputImage, TWaveletFunction>
::GetOutputLowPass() {
  return this->GetOutput(0);
}

template < typename TOutputImage, typename TWaveletFunction>
typename WaveletFrequencyFilterBankForwardGenerator< TOutputImage, TWaveletFunction>
::OutputImagePointer
WaveletFrequencyFilterBankForwardGenerator< TOutputImage, TWaveletFunction>
::GetOutputHighPass() {
  return this->GetOutput(this->m_HighPassSubBands);
}

template < typename TOutputImage, typename TWaveletFunction>
typename WaveletFrequencyFilterBankForwardGenerator< TOutputImage, TWaveletFunction>
::OutputImagePointer
WaveletFrequencyFilterBankForwardGenerator< TOutputImage, TWaveletFunction>
::GetOutputSubBand(unsigned int k) {
  if (k == 0)
    return this->GetOutputLowPass();
  if (k == m_HighPassSubBands)
    return this->GetOutputHighPass();
  return this->GetOutput(k);
}

template < typename TOutputImage, typename TWaveletFunction>
std::vector<typename WaveletFrequencyFilterBankForwardGenerator< TOutputImage, TWaveletFunction>::OutputImagePointer>
WaveletFrequencyFilterBankForwardGenerator< TOutputImage, TWaveletFunction>
::GetOutputs() {
  std::vector<OutputImagePointer> outputList;
  for (unsigned int band = 0; band < this->m_HighPassSubBands + 1; ++band)
    {
    outputList.push_back(this->GetOutputSubBand(band));
    }
  return outputList;
}

template <class TOutputImage, class TWaveletFunction>
void WaveletFrequencyFilterBankForwardGenerator< TOutputImage, TWaveletFunction>
::GenerateData()
// ::ThreadedGenerateData(
//     const OutputImageRegionType& outputRegionForThread,
//     itk::ThreadIdType threadId)
{

  typename TWaveletFunction::Pointer evaluator = TWaveletFunction::New();
  evaluator->SetHighPassSubBands(this->m_HighPassSubBands);

  // Size,Origin and Spacing //
  // Methods inherited from Superclass:GenerateImageSource
  typename OutputImageType::SizeType inputSize = this->GetSize();
  typename OutputImageType::PointType inputOrigin = this->GetOrigin();
  typename OutputImageType::SpacingType inputSpacing = this->GetSpacing();
  typename OutputImageType::SizeType inputHalfSize;
  typename OutputImageType::SpacingType inputSpacingSquare = inputSpacing;
  for (unsigned int i = 0; i < ImageDimension; i++)
    {
    inputSpacingSquare[i] *= inputSpacingSquare[i];
    inputHalfSize[i] = inputSize[i] / 2;
    }

  // Initialize Outputs
  // InputRegionConstIterator inputIt(input, input->GetRequestedRegion());
  std::vector<OutputImagePointer> outputList;
  std::vector<OutputRegionIterator> outputItList;
  for (unsigned int band = 0; band < this->m_HighPassSubBands + 1; ++band)
    {
    outputList.push_back(this->GetOutput(band));
    OutputImagePointer& outputPtr = outputList.back();
    // GenerateImageSource superclass allocates primary output, so use it.(???)
    outputPtr->SetRegions(outputList[0]->GetLargestPossibleRegion());
    outputPtr->Allocate();
    outputPtr->FillBuffer(0);
    outputItList.push_back(OutputRegionIterator(outputPtr,outputPtr->GetRequestedRegion()));
    outputItList.back().GoToBegin();
    }

  // Allocate fake input to iterate
  OutputImagePointer fakeinput = OutputImageType::New();
  fakeinput->SetRegions(outputList.back()->GetLargestPossibleRegion());
  fakeinput->Allocate();
  OutputRegionIterator inputIt(fakeinput, fakeinput->GetRequestedRegion());
  FunctionValueType w2 = 0.0;
  FunctionValueType w = 0.0;
  typename OutputImageType::IndexType index;
  typename OutputImageType::SpacingType w_vector;
  FunctionValueType N_half;
  for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt)
    {
    w2 = 0.0;
    index = inputIt.GetIndex();
    for (unsigned int i = 0; i < ImageDimension; i++)
      {
      if (index[i] <= static_cast<int>(inputHalfSize[i]))
        w_vector[i] = (inputOrigin[i] + inputSpacing[i] * index[i]) /
          static_cast<FunctionValueType>(inputHalfSize[i]);
      else
        w_vector[i] =
          (inputOrigin[i] + inputSpacing[i] * (inputSize[i] - index[i])) /
          static_cast<FunctionValueType>(inputHalfSize[i]);
      w2 += w_vector[i] * w_vector[i];
      }
    w = sqrt(w2);

    itkDebugMacro(<< "w_vector: " << w_vector << " w: " << w
                  << "  inputItIndex: " << inputIt.GetIndex()
                  << " ;; EvaluateHighPassFilter: "
                  << evaluator->EvaluateForwardSubBand(w, m_HighPassSubBands)
                  << " outputIndex: " << outputItList[m_HighPassSubBands].GetIndex());
    // l = 0 is low pass filter, l = m_HighPassSubBands is high-pass filter.
    for (unsigned int l = 0; l < m_HighPassSubBands + 1; ++l)
      {
      outputItList[l].Set(
        outputItList[l].Get() +
        evaluator->EvaluateForwardSubBand(w * this->m_ScaleFactor, l)
        );
      ++outputItList[l];
      }
  }
}
}  // end namespace itk
#endif
