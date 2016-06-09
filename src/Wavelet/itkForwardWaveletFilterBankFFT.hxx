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
template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
ForwardWaveletFilterBankFFT<TInputImage, TOutputImage,
                    TWaveletFunction>::ForwardWaveletFilterBankFFT() {
  this->SetNumberOfRequiredInputs(1);
  this->m_ShrinkFactor = 2;
  this->m_HighPassSubBands = 0;
  this->SetHighPassSubBands(1);
}

template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
void ForwardWaveletFilterBankFFT<TInputImage, TOutputImage,
                    TWaveletFunction>::SetHighPassSubBands(unsigned int k)
{

  if ( m_HighPassSubBands == k )
    {
    return;
    }

  this->m_HighPassSubBands = k;
  this->SetNumberOfRequiredOutputs(k + 1);
  this->Modified();

  for (unsigned int ilevel = 0; ilevel < this->m_HighPassSubBands + 1; ++ilevel)
  {
    this->SetNthOutput(ilevel, this->MakeOutput(ilevel));
  }

}

template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
void ForwardWaveletFilterBankFFT<TInputImage, TOutputImage,
                         TWaveletFunction>::PrintSelf(std::ostream &os,
                                                      Indent indent) const {
  Superclass::PrintSelf(os, indent);

  os << indent << "ShrinkFactor : " << m_ShrinkFactor << std::endl;
}

/*
 * GenerateOutputInformation
 */
template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
void ForwardWaveletFilterBankFFT<TInputImage, TOutputImage,
                         TWaveletFunction>::GenerateOutputInformation() {
  // call the superclass's implementation of this method
  Superclass::GenerateOutputInformation();

  // get pointers to the input and output
  InputImageConstPointer inputPtr = this->GetInput();

  if (!inputPtr) {
    itkExceptionMacro(<< "Input has not been set");
  }

  const typename InputImageType::PointType &inputOrigin = inputPtr->GetOrigin();
  const typename InputImageType::SpacingType &inputSpacing =
      inputPtr->GetSpacing();
  const typename InputImageType::DirectionType &inputDirection =
      inputPtr->GetDirection();
  const typename InputImageType::SizeType &inputSize =
      inputPtr->GetLargestPossibleRegion().GetSize();
  const typename InputImageType::IndexType &inputStartIndex =
      inputPtr->GetLargestPossibleRegion().GetIndex();

  typedef typename OutputImageType::SizeType SizeType;
  typedef typename OutputImageType::IndexType IndexType;

  OutputImagePointer outputPtr;
  typename OutputImageType::PointType outputOrigin;
  typename OutputImageType::SpacingType outputSpacing;
  SizeType outputSize;
  IndexType outputStartIndex;

  // we need to compute the output spacing, the output image size,
  // and the output image start index
  unsigned int low_pass = 0;
  for (unsigned int ilevel = 0; ilevel < this->m_HighPassSubBands + 1; ++ilevel)
    {
    outputPtr = this->GetOutput(ilevel);
    if (!outputPtr) continue;

    for (unsigned int idim = 0; idim < OutputImageType::ImageDimension;
         idim++)
      {
        // Size by half
        outputSize[idim] = static_cast<SizeValueType>(
          std::floor(static_cast<double>(inputSize[idim]) / this->m_ShrinkFactor));
        if (outputSize[idim] < 1) outputSize[idim] = 1;
        // Index by half
        outputStartIndex[idim] = static_cast<IndexValueType>(
          std::ceil(static_cast<double>(inputStartIndex[idim]) / this->m_ShrinkFactor));
        // Spacing is the same!
        outputSpacing[idim] = inputSpacing[idim];
        // Origin
        if (ilevel == low_pass)
          outputOrigin[idim] = inputOrigin[idim];
        else
          outputOrigin[idim] = inputOrigin[idim] +
            outputSpacing[idim] * outputSize[idim] / static_cast<double>(this->m_ShrinkFactor);
      }

    typename OutputImageType::RegionType outputLargestPossibleRegion;
    outputLargestPossibleRegion.SetSize(outputSize);
    outputLargestPossibleRegion.SetIndex(outputStartIndex);

    outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);
    outputPtr->SetOrigin(outputOrigin);
    outputPtr->SetSpacing(outputSpacing);
    outputPtr->SetDirection(inputDirection);  // Output Direction should be same
    // as input.
  }
}

/*
 * GenerateOutputRequestedRegion
 */
template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
void ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>::
    GenerateOutputRequestedRegion(DataObject *refOutput)
{
  // call the superclass's implementation of this method
  Superclass::GenerateOutputRequestedRegion(refOutput);

  // find the index for this output
  unsigned int refLevel =
      static_cast<unsigned int>(refOutput->GetSourceOutputIndex());

  // compute baseIndex and baseSize
  typedef typename OutputImageType::SizeType SizeType;
  typedef typename OutputImageType::IndexType IndexType;
  typedef typename OutputImageType::RegionType RegionType;

  TOutputImage *ptr = itkDynamicCastInDebugMode<TOutputImage *>(refOutput);
  if (!ptr) itkExceptionMacro(<< "Could not cast refOutput to TOutputImage*.");

  unsigned int ilevel, idim;
  unsigned int low_pass = 0;

  if (ptr->GetRequestedRegion() == ptr->GetLargestPossibleRegion()) {
    // set the requested regions for the other outputs to their
    // requested region

  for (unsigned int ilevel = 0; ilevel < this->m_HighPassSubBands + 1; ++ilevel)
      {
      if (ilevel == refLevel) continue;
      if (!this->GetOutput(ilevel)) continue;
      this->GetOutput(ilevel)->SetRequestedRegionToLargestPossibleRegion();
      }
  } else
    {
    // compute requested regions for the other outputs based on
    // the requested region of the reference output
    IndexType outputIndex;
    SizeType outputSize;
    RegionType outputRegion;
    IndexType baseIndex = ptr->GetRequestedRegion().GetIndex();
    SizeType baseSize = ptr->GetRequestedRegion().GetSize();

    for (idim = 0; idim < TOutputImage::ImageDimension; idim++)
      {
      baseIndex[idim] *= static_cast<IndexValueType>(this->m_ShrinkFactor);
      baseSize[idim] *= static_cast<SizeValueType>(this->m_ShrinkFactor);
      }

  for (unsigned int ilevel = 0; ilevel < this->m_HighPassSubBands + 1; ++ilevel)
      {
      if (ilevel == refLevel) continue;
      if (!this->GetOutput(ilevel)) continue;

      for (idim = 0; idim < TOutputImage::ImageDimension; idim++)
        {
        // Index by half.
        outputIndex[idim] = static_cast<IndexValueType>(
          std::ceil(static_cast<double>(baseIndex[idim]) / this->m_ShrinkFactor));
        // Size by half
        outputSize[idim] = static_cast<SizeValueType>(
          std::floor(static_cast<double>(baseSize[idim]) / this->m_ShrinkFactor));
        if (outputSize[idim] < 1) outputSize[idim] = 1;
        outputIndex[idim] = baseIndex[idim];
        }
      }

    outputRegion.SetIndex(outputIndex);
    outputRegion.SetSize(outputSize);

    // make sure the region is within the largest possible region
    outputRegion.Crop(this->GetOutput(ilevel)->GetLargestPossibleRegion());
    // set the requested region
    this->GetOutput(ilevel)->SetRequestedRegion(outputRegion);
  }
}

/**
 * GenerateInputRequestedRegion
 */
template <typename TInputImage, typename TOutputImage,
         typename TWaveletFunction>
void ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>::
GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput());
  if (!inputPtr) {
    itkExceptionMacro(<< "Input has not been set.");
  }

  // compute baseIndex and baseSize
  typedef typename OutputImageType::SizeType SizeType;
  typedef typename OutputImageType::IndexType IndexType;
  typedef typename OutputImageType::RegionType RegionType;

  unsigned int refLevel = 0;
  SizeType baseSize = this->GetOutput(refLevel)->GetRequestedRegion().GetSize();
  IndexType baseIndex =
      this->GetOutput(refLevel)->GetRequestedRegion().GetIndex();
  RegionType baseRegion;

  unsigned int idim;
  for (idim = 0; idim < ImageDimension; idim++) {
    baseIndex[idim] *= static_cast<IndexValueType>(this->m_ShrinkFactor);
    baseSize[idim] *= static_cast<SizeValueType>(this->m_ShrinkFactor);
  }
  baseRegion.SetIndex(baseIndex);
  baseRegion.SetSize(baseSize);

  // make sure the requested region is within the largest possible
  baseRegion.Crop(inputPtr->GetLargestPossibleRegion());

  // set the input requested region
  inputPtr->SetRequestedRegion(baseRegion);
}

template <class TInputImage, class TOutputImage, class TWaveletFunction>
void ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>::
GenerateData()
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
  for (unsigned int i = 0; i < ImageDimension ; i++)
  {
    inputSpacingSquare[i] *= inputSpacingSquare[i];
    inputHalfSize[i] = inputSize[i]/2;
  }

  InputRegionConstIterator inputIt(input, input->GetRequestedRegion());
  std::vector<OutputImagePointer> outputList;
  std::vector<OutputRegionIterator> outItList;
  for (unsigned int ilevel = 0; ilevel < this->m_HighPassSubBands + 1; ++ilevel)
    {
    outputList.push_back(this->GetOutput(ilevel));
    InputImagePointer& outputPtr = outputList.back();
    outputPtr->SetRegions(outputPtr->GetLargestPossibleRegion());
    outputPtr->Allocate();
    outputPtr->FillBuffer(0);
    outItList.push_back(OutputRegionIterator(outputPtr,outputPtr->GetRequestedRegion()));
    }

  OutputImagePointer outputLowPass = outputList[0];
  OutputImagePointer outputHighPass = outputList[1];
  outputLowPass->SetRegions(outputLowPass->GetLargestPossibleRegion());
  outputLowPass->Allocate();
  outputLowPass->FillBuffer(0);

  outputHighPass->SetRegions(outputHighPass->GetLargestPossibleRegion());
  outputHighPass->Allocate();
  outputHighPass->FillBuffer(0);
  OutputRegionIterator outLowIt(outputLowPass, outputLowPass->GetRequestedRegion());
  OutputRegionIterator outHighIt(outputHighPass, outputHighPass->GetRequestedRegion());

  FunctionValueType w2 = 0.0;
  FunctionValueType w = 0.0;
  typename OutputImageType::IndexType index;
  typename OutputImageType::SpacingType w_vector;
  FunctionValueType N_half;
  for( inputIt.GoToBegin(),
       outLowIt.GoToBegin(), outHighIt.GoToBegin();
    !outLowIt.IsAtEnd();
    // SubSample directly here with double ++inputIt
    ++inputIt, ++inputIt, ++outLowIt, ++outHighIt )
    {
    w2 = 0.0;
    index = inputIt.GetIndex();
    // This takes into account origin and spacing information. Sadly, FFT does not update the spacing of the output, so we have to convert evalPoint to frequency.
    for (unsigned int i = 0; i < ImageDimension ; i++)
    {
      if(index[i] <= static_cast<int>(inputHalfSize[i]) )
        w_vector[i] = (inputOrigin[i] + inputSpacing[i] * index[i])/ static_cast<FunctionValueType>(inputHalfSize[i]) ;
      else
        w_vector[i] = (inputOrigin[i] + inputSpacing[i] * (inputSize[i] - index[i]) )/ static_cast<FunctionValueType>(inputHalfSize[i]) ;
      w2 += w_vector[i]*w_vector[i];
    }
    w = sqrt(w2);
    // typename InputImageType::PixelType input_value = inputIt.Get();
    outLowIt.Set( outLowIt.Get() + inputIt.Get() * evaluator->EvaluateForwardLowPassFilter(w));
    outHighIt.Set( outHighIt.Get() + inputIt.Get() * evaluator->EvaluateForwardHighPassFilter(w));

    itkDebugMacro( << "w_vector: " << w_vector << " w2: " << w2 << " ;; Evaluate:: " << evaluator->EvaluateFunction(w2) <<"  inputIt: " << inputIt.GetIndex() <<" outIt: " << outLowIt.GetIndex() );

    }

}
}  // end namespace itk
#endif
