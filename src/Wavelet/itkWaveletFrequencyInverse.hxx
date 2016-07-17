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
#ifndef itkWaveletFrequencyInverse_hxx
#define itkWaveletFrequencyInverse_hxx
#include "itkWaveletFrequencyInverse.h"
#include <itkCastImageFilter.h>
#include "itkImage.h"
#include <algorithm>
#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkExpandImageFilter.h>
#include <itkChangeInformationImageFilter.h>
namespace itk
{
template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
WaveletFrequencyInverse< TInputImage, TOutputImage, TWaveletFilterBank>
::WaveletFrequencyInverse()
  : m_Levels(1),
    m_HighPassSubBands(1),
    m_ScaleFactor(2)
{
  this->SetNumberOfRequiredOutputs(1);
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
std::pair<unsigned int, unsigned int>
WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::InputIndexToLevelBand(unsigned int linear_index)
{
  if (linear_index > this->m_TotalInputs - 1 || linear_index < 0)
    itkExceptionMacro(<< "Failed converting liner index " << linear_index <<
        " to Level,Band pair : out of bounds");
  // Low pass (band = 0).
  if (linear_index == 0 )
    return std::make_pair(this->m_Levels,0);

  unsigned int band = (linear_index - 1) % this->m_HighPassSubBands;
  band = band + 1;
  // note integer division ahead.
  unsigned int level  = (linear_index - band) / this->m_HighPassSubBands + 1;
  return std::make_pair(level, band);
};


template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::SetLevels(unsigned int n)
{
  unsigned int current_inputs = 1 + this->m_Levels * this->m_HighPassSubBands;
  if ( this->m_TotalInputs == current_inputs && this->m_Levels == n )
    {
    return;
    }

  this->m_Levels = n;
  this->m_TotalInputs = 1 + n * this->m_HighPassSubBands;

  this->SetNumberOfRequiredInputs( this->m_TotalInputs );
  this->Modified();

  this->SetNthOutput(0, this->MakeOutput(0));
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::SetHighPassSubBands(unsigned int k)
{
  if ( this->m_HighPassSubBands == k )
    {
    return;
    }
  this->m_HighPassSubBands = k;
  // Trigger setting new number of inputs avoiding code duplication
  this->SetLevels(this->m_Levels);
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void
WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::SetInputs(const std::vector<InputImagePointer> &inputs)
{
  if (inputs.size() != this->m_TotalInputs)
    itkExceptionMacro(<< "Error seting inputs in inverse wavelet. Wrong vector size: " <<
      inputs.size() << " .According to number of levels and bands it should be: " << m_TotalInputs);

  for( unsigned int nin = 0 ; nin < this->m_TotalInputs ; ++nin)
    {
    this->SetNthInput(nin, inputs[nin]);
    }
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void
WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::SetInputLowPass(const InputImagePointer & input_low_pass)
{
  this->SetNthInput(0, input_low_pass);
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void
WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::SetInputsHighPass(const std::vector<InputImagePointer> &inputs)
{
  if (inputs.size() != this->m_TotalInputs - 1)
    itkExceptionMacro(<< "Error seting inputs in inverse wavelet. Wrong vector size: " <<
      inputs.size() << " .According to number of levels and bands it should be: " << m_TotalInputs - 1);

  for( unsigned int nin = 0 ; nin < this->m_TotalInputs - 1 ; ++nin)
    {
    this->SetNthInput(nin + 1, inputs[nin]);
    }

}
template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank >
void WaveletFrequencyInverse< TInputImage, TOutputImage, TWaveletFilterBank>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent <<
    " Levels: " << this->m_Levels <<
    " HighPassSubBands: " << this->m_HighPassSubBands <<
    " TotalInputs: " << this->m_TotalInputs <<
    std::endl;
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
  void WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::GenerateOutputInformation()
{
  // call the superclass's implementation of this method
  Superclass::GenerateOutputInformation();

  // Check  all inputs exist.
  for (unsigned int n_input = 0 ; n_input < this->m_TotalInputs; ++n_input)
    {
    if (!this->GetInput(n_input))
      itkExceptionMacro(<< "Input: " << n_input <<" has not been set");
    }

  // We know inputIndex = 1 has the same size than output. Use it.
  InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput(1));
  const typename InputImageType::PointType &inputOrigin =
    inputPtr->GetOrigin();
  const typename InputImageType::SpacingType &inputSpacing =
    inputPtr->GetSpacing();
  const typename InputImageType::DirectionType &inputDirection =
    inputPtr->GetDirection();
  const typename InputImageType::SizeType &inputSize =
    inputPtr->GetLargestPossibleRegion().GetSize();
  const typename InputImageType::IndexType &inputStartIndex =
    inputPtr->GetLargestPossibleRegion().GetIndex();

  OutputImagePointer outputPtr;
  outputPtr = this->GetOutput(0);

  typename OutputImageType::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize(inputSize);
  outputLargestPossibleRegion.SetIndex(inputStartIndex);

  outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);
  outputPtr->SetOrigin(inputOrigin);
  outputPtr->SetSpacing(inputSpacing);
  outputPtr->SetDirection(inputDirection);
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
  void WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::GenerateOutputRequestedRegion(DataObject *refOutput)
{
  // call the superclass's implementation of this method
  Superclass::GenerateOutputRequestedRegion(refOutput);
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // compute baseIndex and baseSize
  typedef typename OutputImageType::SizeType SizeType;
  typedef typename OutputImageType::IndexType IndexType;
  typedef typename OutputImageType::RegionType RegionType;

  unsigned int refOutput = 0;
  OutputImagePointer outputPtr = this->GetOutput(refOutput);

  if (outputPtr->GetRequestedRegion() == outputPtr->GetLargestPossibleRegion())
    {
    // set the requested regions for inputs to their largest
    for (unsigned int n_input = 0; n_input < this->m_TotalInputs; ++n_input)
      {
      if (!this->GetInput(n_input)) continue;
      InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput(n_input));
      inputPtr->SetRequestedRegionToLargestPossibleRegion();
      }
    }
  else
    {
    IndexType inputIndex;
    SizeType inputSize;
    RegionType inputRegion;
    SizeType baseSize = outputPtr->GetRequestedRegion().GetSize();
    IndexType baseIndex = outputPtr->GetRequestedRegion().GetIndex();
    RegionType baseRegion;
    baseRegion.SetIndex(baseIndex);
    baseRegion.SetSize(baseSize);

    for (unsigned int level = 0; level < this->m_Levels; ++level)
      {
      unsigned int scaleFactorPerLevel = std::pow(this->m_ScaleFactor, level);
      for (unsigned int band = 0 ; band < this->m_HighPassSubBands; ++band)
        {
        unsigned int n_input = 1 + level * this->m_HighPassSubBands + band ;
        if (!this->GetInput(n_input)) continue;
        InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput(n_input));

        for (unsigned int idim = 0; idim < TInputImage::ImageDimension; idim++)
          {
          inputIndex[idim] = baseIndex[idim] * scaleFactorPerLevel;
          inputSize[idim] = baseSize[idim] * scaleFactorPerLevel;
          }

        inputRegion.SetIndex(inputIndex);
        inputRegion.SetSize(inputSize);

        // make sure the region is within the largest possible region
        inputRegion.Crop(inputPtr->GetLargestPossibleRegion());
        // set the requested region
        inputPtr->SetRequestedRegion(inputRegion);

        // Set low pass input
        if(level == this->m_Levels - 1 && band == this->m_HighPassSubBands - 1)
          {
          unsigned int n_input = 0;
          if (!this->GetInput(n_input)) continue;
          InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput(n_input));
          inputRegion.SetIndex(inputIndex);
          inputRegion.SetSize(inputSize);
          inputRegion.Crop(inputPtr->GetLargestPossibleRegion());
          inputPtr->SetRequestedRegion(inputRegion);
          }
        }
      }
    }
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void WaveletFrequencyInverse< TInputImage, TOutputImage, TWaveletFilterBank>
::GenerateData()
{
  this->AllocateOutputs();
  InputImageConstPointer low_pass = this->GetInput(0);

  typedef itk::CastImageFilter<InputImageType, OutputImageType> CastFilterType;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(low_pass);
  castFilter->Update();
  OutputImagePointer low_pass_per_level = castFilter->GetOutput();
  for (int level = this->m_Levels - 1 ; level > -1   ; --level)
    {
    /******* Calculate FilterBank with the right size per level. *****/
    typedef itk::MultiplyImageFilter<OutputImageType> MultiplyFilterType;
    typename WaveletFilterBankType::Pointer filterBank = WaveletFilterBankType::New();
    filterBank->SetHighPassSubBands(this->m_HighPassSubBands);
    filterBank->SetSize(low_pass_per_level->GetLargestPossibleRegion().GetSize() );
    filterBank->Update();

    /** Get new low pass. */
    // TODO: Do I have to upsample before or after?
    // After:::
    typename MultiplyFilterType::Pointer multiplyLowPass = MultiplyFilterType::New();
    multiplyLowPass->SetInput1(filterBank->GetOutputLowPass()) ;
    multiplyLowPass->SetInput2(low_pass_per_level) ;
    multiplyLowPass->InPlaceOn();
    multiplyLowPass->Update();

    typename MultiplyFilterType::Pointer multiplyByConst = MultiplyFilterType::New();
    multiplyByConst->SetInput1(multiplyLowPass->GetOutput());
    multiplyByConst->SetConstant(1 << (ImageDimension + 1)); // 2^dim
    multiplyByConst->Update();
    low_pass_per_level = multiplyByConst->GetOutput();

    /******* set HighPass bands *****/
    std::vector<OutputImagePointer> highPassMasks = filterBank->GetOutputsHighPassBands();
    for(unsigned int band = 0 ; band < this->m_HighPassSubBands ; ++band)
      {
      unsigned int n_input = 1 + level * this->m_HighPassSubBands + band ;

      typename MultiplyFilterType::Pointer multiplyHighBandFilter = MultiplyFilterType::New();
      multiplyHighBandFilter->SetInput1(highPassMasks[band]) ;
      multiplyHighBandFilter->SetInput2(this->GetInput(n_input)) ;
      multiplyHighBandFilter->InPlaceOn();
      multiplyHighBandFilter->Update();

      // Store result summing it to low_pass_per_level
      typedef itk::AddImageFilter<OutputImageType> AddFilterType;
      typename AddFilterType::Pointer addFilter = AddFilterType::New();
      addFilter->SetInput1(low_pass_per_level);
      addFilter->SetInput2(multiplyHighBandFilter->GetOutput());
      addFilter->InPlaceOn();
      addFilter->Update();

      this->UpdateProgress( static_cast< float >( n_input - 1 ) //TODO
        / static_cast< float >( m_TotalInputs ) );
      }

    /******* UpSample for the next Level iteration *****/
    typedef itk::ExpandImageFilter<OutputImageType,OutputImageType> ExpandFilterType;
    typename ExpandFilterType::Pointer expandFilter = ExpandFilterType::New();
    expandFilter->SetInput(low_pass_per_level);
    expandFilter->SetExpandFactors(2);
    expandFilter->Update();
    // low_pass_per_level = expandFilter->GetOutput();
    // Ignore modifications of origin and spacing of expand filters.
    typedef itk::ChangeInformationImageFilter<OutputImageType> ChangeInformationFilterType;
    typename ChangeInformationFilterType::Pointer changeInfoFilter = ChangeInformationFilterType::New();
    changeInfoFilter->SetInput(expandFilter->GetOutput());
    changeInfoFilter->ChangeDirectionOff();
    changeInfoFilter->ChangeRegionOff();
    changeInfoFilter->ChangeSpacingOn();
    changeInfoFilter->ChangeOriginOn();
    changeInfoFilter->UseReferenceImageOn();
    changeInfoFilter->SetReferenceImage(low_pass_per_level.GetPointer()); // Use input image as reference.
    changeInfoFilter->Update();
    low_pass_per_level = changeInfoFilter->GetOutput();

    }
}

} // end namespace itk
#endif
