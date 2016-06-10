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
#include "itkExpandImageFilter.h"
#include "itkChangeInformationImageFilter.h"

namespace itk {
template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
InverseWaveletFilterBankFFT<TInputImage, TOutputImage,
                    TWaveletFunction>::InverseWaveletFilterBankFFT() {
  this->m_ExpandFactor = 2;
  this->m_HighPassSubBands = 0;
  this->SetHighPassSubBands(1);
}

template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage,
                    TWaveletFunction>::SetHighPassSubBands(unsigned int k)
{

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

template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage,
                         TWaveletFunction>::PrintSelf(std::ostream &os,
                                                      Indent indent) const {
  Superclass::PrintSelf(os, indent);

  os << indent << "ExpandFactor : "    << this->m_ExpandFactor <<
               " ; HighPassSubBands : " << this->m_HighPassSubBands << std::endl;
}

/* ******* Set Input *****/
template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::SetInputLowPass(InputImagePointer imgP) {
  return this->SetNthInput(0, imgP);
}

template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::SetInputHighPass(InputImagePointer imgP) {
  return this->SetNthInput(this->m_HighPassSubBands, imgP);
}

template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::SetInputSubBand(unsigned int k, InputImagePointer imgP) {
  if (k==0) return this->SetInputLowPass(imgP);
  if (k==m_HighPassSubBands) return this->SetInputHighPass(imgP);
  return this->SetNthInput(k, imgP);
}

template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
void
InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::SetInputs(std::vector<InputImagePointer> inputList) {
  for (unsigned int ilevel = 0; ilevel < this->m_HighPassSubBands + 1; ++ilevel)
    {
    this->SetInputSubBand(ilevel, inputList[ilevel]);
    }
}
/*
 * GenerateOutputInformation
 */
template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage,
                         TWaveletFunction>::GenerateOutputInformation() {
  // call the superclass's implementation of this method
  Superclass::GenerateOutputInformation();

  // get pointers to the low pass(just choose one) input and output.
  InputImageConstPointer inputPtr = this->GetInput(0);

  if (!inputPtr)
    itkExceptionMacro(<< "Input has not been set");

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

  outputPtr = this->GetOutput();

  for (unsigned int idim = 0; idim < OutputImageType::ImageDimension;
    idim++)
    {
    // Size double
    outputSize[idim] = static_cast<SizeValueType>(
      std::floor(static_cast<double>(inputSize[idim]) * this->m_ExpandFactor));
    // Index double
    outputStartIndex[idim] = static_cast<IndexValueType>(
      std::ceil(static_cast<double>(inputStartIndex[idim]) * this->m_ExpandFactor));
    // Spacing is the same!
    outputSpacing[idim] = inputSpacing[idim];
    // Origin
    outputOrigin[idim] = inputOrigin[idim];
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

/**
 * GenerateInputRequestedRegion
 */
template <typename TInputImage, typename TOutputImage,
         typename TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>::
GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput(0));
  if (!inputPtr)
    itkExceptionMacro(<< "Input has not been set.");

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
    baseSize[idim] = static_cast<SizeValueType>(
      std::floor(static_cast<double>(baseSize[idim]) / this->m_ExpandFactor));
    baseIndex[idim] = static_cast<IndexValueType>(
      std::ceil(static_cast<double>(baseIndex[idim]) / this->m_ExpandFactor));
  }
  baseRegion.SetIndex(baseIndex);
  baseRegion.SetSize(baseSize);

  // make sure the requested region is within the largest possible
  baseRegion.Crop(inputPtr->GetLargestPossibleRegion());

  // set the input requested region
  inputPtr->SetRequestedRegion(baseRegion);
}

template <class TInputImage, class TOutputImage, class TWaveletFunction>
void InverseWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>::
GenerateData()
// ::ThreadedGenerateData(
//     const OutputImageRegionType& outputRegionForThread,
//     itk::ThreadIdType threadId)
{
  std::vector<InputImageConstPointer> inputs;
  for (unsigned int l = 0; l < m_HighPassSubBands + 1 ; ++l)
    {
    inputs.push_back(this->GetInput(l));
    }
  InputImageConstPointer input = inputs[0];

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
  std::vector<OutputImagePointer> preDownSampledList;
  std::vector<OutputRegionIterator> preDownSampledItList;
  std::vector<OutputImagePointer> outputList;
  std::vector<OutputRegionIterator> outItList;
  for (unsigned int ilevel = 0; ilevel < this->m_HighPassSubBands + 1; ++ilevel)
    {
    preDownSampledList.push_back(InputImageType::New());
    InputImagePointer& preDownSampledPtr = preDownSampledList.back();
    preDownSampledPtr->SetRegions(input->GetLargestPossibleRegion());
    preDownSampledPtr->Allocate();
    preDownSampledPtr->FillBuffer(0);
    preDownSampledItList.push_back(OutputRegionIterator(preDownSampledPtr,
        preDownSampledPtr->GetRequestedRegion()));
    preDownSampledItList.back().GoToBegin();

    outputList.push_back(this->GetOutput(ilevel));
    InputImagePointer& outputPtr = outputList.back();
    outputPtr->SetRegions(outputPtr->GetLargestPossibleRegion());
    outputPtr->Allocate();
    outputPtr->FillBuffer(0);
    outItList.push_back(OutputRegionIterator(outputPtr,outputPtr->GetRequestedRegion()));
    outItList.back().GoToBegin();
    }

  FunctionValueType w2 = 0.0;
  FunctionValueType w = 0.0;
  typename OutputImageType::IndexType index;
  typename OutputImageType::SpacingType w_vector;
  FunctionValueType N_half;
  for( inputIt.GoToBegin();
    !inputIt.IsAtEnd();
    ++inputIt )
    {
    w2 = 0.0;
    index = inputIt.GetIndex();
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
      "  inputItIndex: " << inputIt.GetIndex() <<
      " ;; l0:EvaluateLowPassFilter: " << evaluator->EvaluateInverseSubBand(w,0)  <<
      " preDownSampledIndex: " << preDownSampledItList[0].GetIndex() );
    // l = 0 is low pass filter, l = m_HighPassSubBands is high-pass filter.
    for (unsigned int l = 0; l < m_HighPassSubBands + 1 ; ++l)
      {
      preDownSampledItList[l].Set( preDownSampledItList[l].Get() +
        inputIt.Get() * evaluator->EvaluateInverseSubBand(w,l));
      ++preDownSampledItList[l];
      }
    }

  // Downsample:
  // Not really interested in the output information of the shrinker,
  // so we change it to that of the output which was set in GenerateOutputInformation.
  typedef itk::ExpandImageFilter<OutputImageType,OutputImageType> ExpanderType;
  unsigned int factors[ImageDimension];
  for (unsigned int i = 0 ; i< ImageDimension ; ++i)
    factors[i] = this->m_ExpandFactor;

  typedef itk::ChangeInformationImageFilter<OutputImageType> ChangeInformationType;

  for (unsigned int ilevel = 0; ilevel < this->m_HighPassSubBands + 1; ++ilevel)
    {
    typename ExpanderType::Pointer shrinker = ExpanderType::New();
    shrinker->SetExpandFactors(factors);
    shrinker->SetInput(preDownSampledList[ilevel]);
    typename ChangeInformationType::Pointer changerInfo = ChangeInformationType::New();
    changerInfo->ChangeAll();
    // This actually copy the information set in GenerateOutputInformation().
    changerInfo->SetReferenceImage(outputList[ilevel]);
    changerInfo->SetInput(shrinker->GetOutput());
    changerInfo->Update();
    this->GraftNthOutput(ilevel, changerInfo->GetOutput());
    }
}
}  // end namespace itk
#endif
