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
#ifndef itkForwardWaveletFFT_hxx
#define itkForwardWaveletFFT_hxx
#include "itkForwardWaveletFFT.h"
#include <itkCastImageFilter.h>
#include "visualize_functions.h" //TODO delete
#include "itkInverseFFTImageFilter.h" // TODO delete
#include "itkImage.h"
#include <algorithm>
#include <itkMultiplyImageFilter.h>
#define __my_debug__ 1
namespace itk
{
template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
ForwardWaveletFFT< TInputImage, TOutputImage, TWaveletFilterBank>
::ForwardWaveletFFT()
  : m_Levels(1),
    m_HighPassSubBands(1),
    m_TotalOutputs(1)
{
  this->SetNumberOfRequiredInputs(1);
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
unsigned int ForwardWaveletFFT<TInputImage, TOutputImage, TWaveletFilterBank>
::ComputeMaxNumberOfLevels(typename InputImageType::SizeType& input_size)
{
  std::array<unsigned int, ImageDimension> exponent_per_axis;
  for (unsigned int axis = 0 ; axis < ImageDimension ; ++axis)
  {
    size_t size_axis = input_size[axis];
    if (size_axis < 2)
    {
      exponent_per_axis[axis] = 1;
      continue;
    }
    double exponent = std::log(size_axis) / std::log(2.0) ;
    // check that exponent is integer: the fractional part is 0
    double int_part;
    if (std::modf(exponent, &int_part ) == 0.0 )
    {
      exponent_per_axis[axis] = static_cast<unsigned int>(exponent);
      continue;
    }
    else
    {
      exponent_per_axis[axis] = 1;
      continue;
    }
  }
  // return the min_element of array (1 if any size is not power of 2)
  return *std::min_element(exponent_per_axis.begin(), exponent_per_axis.end());
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void ForwardWaveletFFT<TInputImage, TOutputImage, TWaveletFilterBank>
::SetLevels(unsigned int n)
{
  // outputs = lowpass:1 + high_pass:levels * bands
  unsigned int current_outputs = 1 + this->m_Levels * this->m_HighPassSubBands;
  if ( this->m_TotalOutputs == current_outputs && this->m_Levels == n )
    {
    return;
    }

  this->m_Levels = n;
  this->m_TotalOutputs = 1 + n * this->m_HighPassSubBands;

  this->SetNumberOfRequiredOutputs( this->m_TotalOutputs );
  this->Modified();

  for (unsigned int n_output = 0; n_output < this->m_TotalOutputs; ++n_output)
  {
    this->SetNthOutput(n_output, this->MakeOutput(n_output));
  }
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void ForwardWaveletFFT<TInputImage, TOutputImage, TWaveletFilterBank>
::SetHighPassSubBands(unsigned int k)
{
  if ( this->m_HighPassSubBands == k )
    {
    return;
    }
  this->m_HighPassSubBands = k;
  // Trigger setting new number of outputs avoiding code duplication
  this->SetLevels(this->m_Levels);
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank >
void ForwardWaveletFFT< TInputImage, TOutputImage, TWaveletFilterBank>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent <<
    " Levels: " << this->m_Levels <<
    " HighPassSubBands: " << this->m_HighPassSubBands <<
    " TotalOutputs: " << this->m_TotalOutputs <<
    std::endl;
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void ForwardWaveletFFT< TInputImage, TOutputImage, TWaveletFilterBank>
::GenerateData()
{
  InputImageConstPointer input = this->GetInput();
  // this->AllocateOutputs();

  std::vector<OutputImagePointer> outputList;
  std::vector<OutputRegionIterator> outputItList;
  // typedef itk::CastImageFilter<InputImageType, OutputImageType> CastFilterType;
  // typename CastFilterType::Pointer castFilter = CastFilterType::New();
  // castFilter->SetInput(input);
  // castFilter->Update();
  // First level:
  unsigned int ilevel = 0;
  typename WaveletFilterBankType::Pointer filterBank = WaveletFilterBankType::New();
  filterBank->SetHighPassSubBands(this->m_HighPassSubBands);
  filterBank->SetSize(input->GetLargestPossibleRegion().GetSize() );
  filterBank->Update();
  typedef itk::MultiplyImageFilter<OutputImageType> MultiplyFilterType;
  typename MultiplyFilterType::Pointer multiplyLowFilter = MultiplyFilterType::New();
  /* The first output is the low_pass image.
   * It returns the input image multiplied by the low-pass filter, recursively at each level.
   */
  OutputImagePointer outputLowPass = this->GetOutput(0);
  outputLowPass->SetRegions( outputLowPass->GetRequestedRegion() );
  outputLowPass->Allocate();

  multiplyLowFilter->SetInput1(filterBank->GetOutputLowPass()) ;
  multiplyLowFilter->SetInput2(input) ;
  multiplyLowFilter->InPlaceOn();
  multiplyLowFilter->GraftOutput(outputLowPass);
  // multiplyLowFilter->Modified();
  // multiplyLowFilter->UpdateLargestPossibleRegion();
  multiplyLowFilter->Update();
  this->GraftNthOutput(0, multiplyLowFilter->GetOutput());

  /* The other outputs are the high-pass at each level.( and each sub-band)
   * Multiply HighPass (bands too) waveletFilter bank with input image,
   * and save the result.
   */
  std::vector<OutputImagePointer> highPassImages = filterBank->GetOutputsHighPassBands();
  std::vector<OutputImagePointer> outputMultiplyHigh;
  for(unsigned int band = 0 ; band < this->m_HighPassSubBands ; ++band)
  {
    unsigned int NOutput = 1 + ilevel + band;
    OutputImagePointer outputHighSubBand = this->GetOutput(NOutput);
    outputHighSubBand->SetRegions( outputHighSubBand->GetLargestPossibleRegion() );
    outputHighSubBand->Allocate();

    typename MultiplyFilterType::Pointer multiplyHighBandFilter = MultiplyFilterType::New();
    multiplyHighBandFilter->SetInput1(highPassImages[band]) ;
    multiplyHighBandFilter->SetInput2(input) ;
    multiplyHighBandFilter->InPlaceOn();
    multiplyHighBandFilter->GraftOutput(outputHighSubBand);
    // multiplyHighBandFilter->Modified();
    multiplyHighBandFilter->Update();
    this->GraftNthOutput(NOutput, multiplyHighBandFilter->GetOutput());
  }
  ++ilevel;

  // std::cout << outputLowPass->GetLargestPossibleRegion() << std::endl;
  // for (unsigned int level = 0; level < this->m_Levels; ++level)
  // {
  //   // Low pass is the first.
  //   if(level == 0)
  //   {
  //     this->GraftNthOutput( 0, outputLowPass );
  //   }
  //
    // Calculate the index of the first output of this level.
    // +1 because we only return one (the first index) low-pass image in the whole pyramid.
  //   unsigned int first_output_current_level =
  //     static_cast<unsigned int>( 1 + this->m_HighPassSubBands * level );
  //
  //   // Bands per level
  //   for (unsigned int band = 0 ; band < this->m_HighPassSubBands; ++band)
  //   {
  //     this->GraftNthOutput( first_output_current_level + band , outputsHighPerLevel[level][band] );
  //   }
  // }
//   //TODO: 
//   for (; ilevel < this->m_Levels; ++ilevel)
//     {
//     for(unsigned int band = 0 ; band < this->m_HighPassSubBands + 1 ; ++band)
//       {
//       typename WaveletFilterBankType::Pointer filterBank = WaveletFilterBankType::New();
//       filterBank->SetHighPassSubBands(this->m_HighPassSubBands);
//       filterBank->SetInput(outputsHighPerLevel[ilevel - 1][band]);
//       filterBank->Update();
//       std::vector<OutputImagePointer> outputsPerBand = filterBank->GetOutputs();
//       if (band == 0)
//         outputsHighPerLevel.push_back(outputsPerBand);
//       else
//       outputsHighPerLevel[ilevel].insert(
//         outputsHighPerLevel[ilevel].end(),
//         outputsPerBand.begin(),
//         outputsPerBand.end());
//       }
//     }
//
//   for (unsigned int level = 0, n_out_old_level = 0; level < this->m_Levels; ++level)
//     {
//     unsigned int n_out_current_level = static_cast<unsigned int>(std::pow(this->m_HighPassSubBands + 1, level + 1));
//     for (unsigned int i = 0 ; i < n_out_current_level; ++i)
//       {
//       this->GraftNthOutput(n_out_old_level + i, outputsHighPerLevel[level][i] );
// #if __my_debug__ // TODO DELETE
//       std::cout << outputsHighPerLevel[level][i]->GetBufferedRegion() << std::endl;
//       typedef itk::InverseFFTImageFilter<OutputImageType, itk::Image<double, ImageDimension> > InverseFFTFilterType;
//       typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();
//       inverseFFT->SetInput(outputsHighPerLevel[level][i]);
//       inverseFFT->Update();
//       std::cout <<"NBands:" << this->m_HighPassSubBands + 1 <<  " ;; (Level, output) = " << level << " , " << i << " Sublevel: " << i/(this->m_HighPassSubBands + 1) << " SubOutput: " << i%(this->m_HighPassSubBands +1)<< std::endl;
//       visualize::VisualizeITKImage(inverseFFT->GetOutput());
// #endif
//       }
// #if __my_debug__ // TODO DELETE
//   std::cout <<"m_Levels : " << this->m_Levels << " Level: " << level << " size: " << outputsHighPerLevel[level].size() <<std::endl;
// #endif
//     n_out_old_level = n_out_current_level;
//     }

}

} // end namespace itk
#endif
