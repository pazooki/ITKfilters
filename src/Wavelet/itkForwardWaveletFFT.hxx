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
#define __my_debug__ 1
namespace itk
{
template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
ForwardWaveletFFT< TInputImage, TOutputImage, TWaveletFilterBank>
::ForwardWaveletFFT()
{
  this->SetNumberOfRequiredInputs(1);
  this->m_Levels = 1;
  this->m_HighPassSubBands = 1;
  this->m_TotalOutputs = 1;
}

template <typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void ForwardWaveletFFT<TInputImage, TOutputImage, TWaveletFilterBank>
::SetLevels(unsigned int n)
{
  unsigned int current_outputs = static_cast<unsigned int>(
     (1 - static_cast<int>(std::pow( this->m_HighPassSubBands + 1, this->m_Levels + 1 )) ) / (1 - static_cast<int>(m_HighPassSubBands + 1)) - 1 );

  if ( this->m_TotalOutputs == current_outputs && this->m_Levels == n )
    {
    return;
    }

  this->m_Levels = n;
  this->m_TotalOutputs = static_cast<unsigned int>(
     (1 - static_cast<int>(std::pow( this->m_HighPassSubBands + 1, this->m_Levels + 1 )) ) / (1 - static_cast<int>(m_HighPassSubBands + 1)) - 1 );

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
  this->AllocateOutputs();

  std::vector<std::vector<OutputImagePointer>> outputsPerLevel;
  std::vector<OutputImagePointer> outputList;
  std::vector<OutputRegionIterator> outputItList;
  typedef itk::CastImageFilter<InputImageType, OutputImageType> CastFilterType;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(input);
  // First level:
  unsigned int ilevel = 0;
  typename WaveletFilterBankType::Pointer filterBank = WaveletFilterBankType::New();
  filterBank->SetHighPassSubBands(this->m_HighPassSubBands);
  filterBank->SetInput(castFilter->GetOutput());
  filterBank->Update();
  outputsPerLevel.push_back( filterBank->GetOutputs() );
  ++ilevel;

  for (; ilevel < this->m_Levels; ++ilevel)
    {
    for(unsigned int band = 0 ; band < this->m_HighPassSubBands + 1 ; ++band)
      {
      typename WaveletFilterBankType::Pointer filterBank = WaveletFilterBankType::New();
      filterBank->SetHighPassSubBands(this->m_HighPassSubBands);
      filterBank->SetInput(outputsPerLevel[ilevel - 1][band]);
      filterBank->Update();
      std::vector<OutputImagePointer> outputsPerBand = filterBank->GetOutputs();
      if (band == 0)
        outputsPerLevel.push_back(outputsPerBand);
      else
      outputsPerLevel[ilevel].insert(
        outputsPerLevel[ilevel].end(),
        outputsPerBand.begin(),
        outputsPerBand.end());
      }
    }

  for (unsigned int level = 0, n_out_old_level = 0; level < this->m_Levels; ++level)
    {
    unsigned int n_out_current_level = static_cast<unsigned int>(std::pow(this->m_HighPassSubBands + 1, level + 1));
    for (unsigned int i = 0 ; i < n_out_current_level; ++i)
      {
      this->GraftNthOutput(n_out_old_level + i, outputsPerLevel[level][i] );
#if __my_debug__ // TODO DELETE
      std::cout << outputsPerLevel[level][i]->GetBufferedRegion() << std::endl;
      typedef itk::InverseFFTImageFilter<OutputImageType, itk::Image<double, ImageDimension> > InverseFFTFilterType;
      typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();
      inverseFFT->SetInput(outputsPerLevel[level][i]);
      inverseFFT->Update();
      std::cout <<"NBands:" << this->m_HighPassSubBands + 1 <<  " ;; (Level, output) = " << level << " , " << i << " Sublevel: " << i/(this->m_HighPassSubBands + 1) << " SubOutput: " << i%(this->m_HighPassSubBands +1)<< std::endl;
      visualize::VisualizeITKImage(inverseFFT->GetOutput());
#endif
      }
#if __my_debug__ // TODO DELETE
  std::cout <<"m_Levels : " << this->m_Levels << " Level: " << level << " size: " << outputsPerLevel[level].size() <<std::endl;
#endif
    n_out_old_level = n_out_current_level;
    }


}

} // end namespace itk
#endif
