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
#ifndef itkInverseWaveletFilterBankFFT_hxx
#define itkInverseWaveletFilterBankFFT_hxx
#include "itkInverseWaveletFilterBankFFT.h"
#include "itkNumericTraits.h"

namespace itk {
template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::InverseWaveletFilterBankFFT() {
  this->m_ExpandFactor = 2;
  this->m_HighPassSubBands = 0;
  this->SetHighPassSubBands(1);
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::SetHighPassSubBands(unsigned int k) {

  this->SetNumberOfRequiredOutputs(1);
  if ( m_HighPassSubBands == k )
    {
    return;
    }

  this->m_HighPassSubBands = k;
  this->SetNumberOfRequiredInputs(k + 1);
  this->Modified();

  this->SetNthOutput(0, this->MakeOutput(0));

}

template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::PrintSelf(std::ostream &os, Indent indent) const {
  Superclass::PrintSelf(os, indent);

  os << indent << "ExpandFactor : "    << this->m_ExpandFactor <<
               " ; HighPassSubBands : " << this->m_HighPassSubBands << std::endl;
}

/* ******* Set Input *****/
template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::SetInputLowPass(InputImagePointer imgP) {
  return this->SetNthInput(0, imgP);
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::SetInputHighPass(InputImagePointer imgP) {
  return this->SetNthInput(this->m_HighPassSubBands, imgP);
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::SetInputSubBand(unsigned int k, InputImagePointer imgP) {
  if (k==0) return this->SetInputLowPass(imgP);
  if (k==m_HighPassSubBands) return this->SetInputHighPass(imgP);
  return this->SetNthInput(k, imgP);
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::SetInputs(std::vector<InputImagePointer> inputList) {
  for (unsigned int ilevel = 0; ilevel < this->m_HighPassSubBands + 1; ++ilevel)
    {
    this->SetInputSubBand(ilevel, inputList[ilevel]);
    }
}

template <class TInputImage, class TOutputImage, class TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::GenerateData()
// ::ThreadedGenerateData(
//     const OutputImageRegionType& outputRegionForThread,
//     itk::ThreadIdType threadId)
{
  this->AllocateOutputs();
  InputImageConstPointer input = this->GetInput(0);
  // Set Sizes (input and output has the same size)
  typename OutputImageType::SizeType inputSize = input->GetLargestPossibleRegion().GetSize();
  typename OutputImageType::PointType inputOrigin = input->GetOrigin();
  typename OutputImageType::SpacingType inputSpacing = input->GetSpacing();
  typename InputImageType::SizeType inputHalfSize;
  typename InputImageType::SpacingType inputSpacingSquare = input->GetSpacing();
  for (unsigned int i = 0; i < ImageDimension ; i++)
  {
    inputSpacingSquare[i] *= inputSpacingSquare[i];
    inputHalfSize[i] = inputSize[i]/2;
  }
  std::vector<InputImageConstPointer> inputs;
  std::vector<InputRegionConstIterator> inputsIt;
  for (unsigned int l = 0; l < m_HighPassSubBands + 1 ; ++l)
    {
    inputs.push_back(this->GetInput(l));
    inputsIt.push_back(InputRegionConstIterator(inputs.back(),inputs.back()->GetRequestedRegion()));
    inputsIt.back().GoToBegin();
    }

  typename TWaveletFunction::Pointer evaluator = TWaveletFunction::New();
  evaluator->SetHighPassSubBands(this->m_HighPassSubBands);

  OutputImagePointer outputPtr = this->GetOutput();
  OutputRegionIterator outputIt(outputPtr, outputPtr->GetRequestedRegion());

  FunctionValueType w2 = 0.0;
  FunctionValueType w = 0.0;
  typename OutputImageType::IndexType index;
  typename OutputImageType::SpacingType w_vector;
  FunctionValueType N_half;
  for( outputIt.GoToBegin();
    !outputIt.IsAtEnd();
    ++outputIt )
    {
    w2 = 0.0;
    index = outputIt.GetIndex();
    for (unsigned int i = 0; i < ImageDimension ; i++)
    {
      if(index[i] <= static_cast<int>(inputHalfSize[i]) )
        w_vector[i] = (inputOrigin[i] + inputSpacing[i] * index[i]) /
          static_cast<FunctionValueType>(inputHalfSize[i]) ;
      else
        w_vector[i] = (inputOrigin[i] + inputSpacing[i] * (inputSize[i] - index[i])) /
          static_cast<FunctionValueType>(inputHalfSize[i]) ;
      w2 += w_vector[i]*w_vector[i];
    }
    w = sqrt(w2);

    itkDebugMacro( << "w_vector: " << w_vector << " w: " << w <<
      " ;; l0:EvaluateLowPassFilter: " << evaluator->EvaluateInverseSubBand(w,0)  <<
      " outputIndex: " << outputIt.GetIndex() );
    // l = 0 is low pass filter, l = m_HighPassSubBands is high-pass filter.
    for (unsigned int l = 0; l < m_HighPassSubBands + 1 ; ++l)
      {
      outputIt.Set( outputIt.Get() + inputsIt[l].Get() * evaluator->EvaluateInverseSubBand(w,l));
      ++inputsIt[l];
      }
    }

}
}  // end namespace itk
#endif
