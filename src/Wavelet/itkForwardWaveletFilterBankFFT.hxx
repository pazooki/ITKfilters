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
#include "itkShrinkImageFilter.h"
#include "itkChangeInformationImageFilter.h"

namespace itk {
template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
ForwardWaveletFilterBankFFT<TInputImage, TOutputImage,
                    TWaveletFunction>::ForwardWaveletFilterBankFFT() {
  this->m_ShrinkFactor = 2;
  this->m_HighPassSubBands = 0;
  this->SetHighPassSubBands(1);
}

template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
void ForwardWaveletFilterBankFFT<TInputImage, TOutputImage,
                    TWaveletFunction>::SetHighPassSubBands(unsigned int k)
{

  this->SetNumberOfRequiredInputs(1);
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

  os << indent << "ShrinkFactor : "    << this->m_ShrinkFactor <<
               " ; HighPassSubBands : " << this->m_HighPassSubBands << std::endl;
}

/* ******* Get Outputs *****/
template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
typename ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::OutputImagePointer
ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::GetOutputLowPass() {
  return this->GetOutput(0);
}
template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
typename ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::OutputImagePointer
ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::GetOutputHighPass() {
  return this->GetOutput(this->m_HighPassSubBands);
}
template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
typename ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::OutputImagePointer
ForwardWaveletFilterBankFFT<TInputImage, TOutputImage, TWaveletFunction>
::GetOutputSubBand(unsigned int k) {
  if (k==0) return this->GetOutputLowPass();
  if (k==m_HighPassSubBands) return this->GetOutputHighPass();
  return this->GetOutput(k);
}
template <typename TInputImage, typename TOutputImage,
          typename TWaveletFunction>
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
        outputOrigin[idim] = inputOrigin[idim];
        // if (ilevel == low_pass)
        //   outputOrigin[idim] = inputOrigin[idim];
        // else
        //   outputOrigin[idim] = inputOrigin[idim] +
        //     outputSpacing[idim] * outputSize[idim] / static_cast<double>(this->m_ShrinkFactor);
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

      outputRegion.SetIndex(outputIndex);
      outputRegion.SetSize(outputSize);

      // make sure the region is within the largest possible region
      outputRegion.Crop(this->GetOutput(ilevel)->GetLargestPossibleRegion());
      // set the requested region
      this->GetOutput(ilevel)->SetRequestedRegion(outputRegion);
      }
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
      " ;; l0:EvaluateLowPassFilter: " << evaluator->EvaluateForwardSubBand(w,0)  <<
      " preDownSampledIndex: " << preDownSampledItList[0].GetIndex() );
    // l = 0 is low pass filter, l = m_HighPassSubBands is high-pass filter.
    for (unsigned int l = 0; l < m_HighPassSubBands + 1 ; ++l)
      {
      preDownSampledItList[l].Set( preDownSampledItList[l].Get() +
        inputIt.Get() * evaluator->EvaluateForwardSubBand(w,l));
      ++preDownSampledItList[l];
      }
    }

  // Downsample:
  // Not really interested in the output information of the shrinker,
  // so we change it to that of the output which was set in GenerateOutputInformation.
  typedef itk::ShrinkImageFilter<OutputImageType,OutputImageType> ShrinkerType;
  unsigned int factors[ImageDimension];
  for (unsigned int i = 0 ; i< ImageDimension ; ++i)
    factors[i] = this->m_ShrinkFactor;

  typedef itk::ChangeInformationImageFilter<OutputImageType> ChangeInformationType;

  for (unsigned int ilevel = 0; ilevel < this->m_HighPassSubBands + 1; ++ilevel)
    {
    typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
    shrinker->SetShrinkFactors(factors);
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
