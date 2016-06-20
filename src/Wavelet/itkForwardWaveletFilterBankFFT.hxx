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
#ifndef itkForwardWaveletFilterBankFFT_hxx
#define itkForwardWaveletFilterBankFFT_hxx
#include "itkForwardWaveletFilterBankFFT.h"
#include "itkNumericTraits.h"

namespace itk {
template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::ForwardWaveletFilterBankFFT() {
  this->m_ShrinkFactor = 2;
  this->m_HighPassSubBands = 0;
  this->m_ScaleFactor = 1;
  this->SetHighPassSubBands(1);
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
void ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::SetHighPassSubBands(unsigned int k)
{

  this->SetNumberOfRequiredInputs(1);
  if ( m_HighPassSubBands == k ) return;

  this->m_HighPassSubBands = k;
  this->SetNumberOfRequiredOutputs(k + 1);
  this->Modified();

  for (unsigned int ilevel = 0; ilevel < this->m_HighPassSubBands + 1; ++ilevel)
    this->SetNthOutput(ilevel, this->MakeOutput(ilevel));

}

template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
void ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::PrintSelf(std::ostream &os, Indent indent) const {
  Superclass::PrintSelf(os, indent);

  os << indent << "ShrinkFactor : " << this->m_ShrinkFactor
     << " ; HighPassSubBands : " << this->m_HighPassSubBands << std::endl;
}

/* ******* Get Outputs *****/
template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
typename ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::OutputImagePointer
ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::GetOutputLowPass() {
  return this->GetOutput(0);
}
template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
typename ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::OutputImagePointer
ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::GetOutputHighPass() {
  return this->GetOutput(this->m_HighPassSubBands);
}
template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
typename ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::OutputImagePointer
ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::GetOutputSubBand(unsigned int k) {
  if (k == 0)
    return this->GetOutputLowPass();
  if (k == m_HighPassSubBands)
    return this->GetOutputHighPass();
  return this->GetOutput(k);
}
template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
std::vector<typename ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>::OutputImagePointer>
ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::GetOutputs() {
  std::vector<OutputImagePointer> outputList;
  for (unsigned int ilevel = 0; ilevel < this->m_HighPassSubBands + 1; ++ilevel)
    {
    outputList.push_back(this->GetOutputSubBand(ilevel));
    }
  return outputList;
}

template <class TInputImage, class TOutputImage, class TWaveletFunction>
void ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::GenerateData()
// ::ThreadedGenerateData(
//     const OutputImageRegionType& outputRegionForThread,
//     itk::ThreadIdType threadId)
{
  InputImageConstPointer input = this->GetInput();

  typename TWaveletFunction::Pointer evaluator = TWaveletFunction::New();
  evaluator->SetHighPassSubBands(this->m_HighPassSubBands);

  typename OutputImageType::SizeType inputSize = input->GetLargestPossibleRegion().GetSize();
  typename OutputImageType::PointType inputOrigin = input->GetOrigin();
  typename OutputImageType::SpacingType inputSpacing = input->GetSpacing();
  typename InputImageType::SizeType inputHalfSize;
  typename InputImageType::SpacingType inputSpacingSquare = input->GetSpacing();
  for (unsigned int i = 0; i < ImageDimension; i++)
    {
    inputSpacingSquare[i] *= inputSpacingSquare[i];
    inputHalfSize[i] = inputSize[i] / 2;
    }

  InputRegionConstIterator inputIt(input, input->GetRequestedRegion());
  std::vector<OutputImagePointer> outputList;
  std::vector<OutputRegionIterator> outputItList;
  for (unsigned int ilevel = 0; ilevel < this->m_HighPassSubBands + 1; ++ilevel)
    {
    outputList.push_back(this->GetOutput(ilevel));
    InputImagePointer& outputPtr = outputList.back();
    outputPtr->SetRegions(outputPtr->GetLargestPossibleRegion());
    outputPtr->Allocate();
    outputPtr->FillBuffer(0);
    outputItList.push_back(OutputRegionIterator(outputPtr,outputPtr->GetRequestedRegion()));
    outputItList.back().GoToBegin();
    }

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
                  << " ;; l0:EvaluateLowPassFilter: "
                  << evaluator->EvaluateForwardSubBand(w, 0)
                  << " outputIndex: " << outputItList[0].GetIndex());
    // l = 0 is low pass filter, l = m_HighPassSubBands is high-pass filter.
    for (unsigned int l = 0; l < m_HighPassSubBands + 1; ++l)
      {
      outputItList[l].Set(
        outputItList[l].Get() +
        inputIt.Get() *
          evaluator->EvaluateForwardSubBand(w * this->m_ScaleFactor, l));
      ++outputItList[l];
      }
  }
}
}  // end namespace itk
#endif
