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
#ifndef itkVowWavelet_hxx
#define itkVowWavelet_hxx

#include <cmath>
#include "itkMath.h"
#include "itkVowWavelet.h"

namespace itk
{
template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
VowWavelet< TFunctionValue, VImageDimension, TInput >
::VowWavelet()
 :m_Kappa(0.75)
{}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
VowWavelet< TFunctionValue, VImageDimension, TInput >
::~VowWavelet()
{}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename VowWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
VowWavelet< TFunctionValue, VImageDimension, TInput >
::Evaluate(const TInput & position) const
{
  TFunctionValue freq_in_hz = 0.0;
  for (unsigned int i = 0; i < VImageDimension ; i++)
    freq_in_hz += position[i]*position[i];

  return this->EvaluateFunction(freq_in_hz);
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename VowWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
VowWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateFunction(const FunctionValueType & freq_in_hz) const
{
  // freq_in_rad_per_sec = freq_in_hz * 2 * pi
  if( freq_in_hz >= 1/8.0 && freq_in_hz < 1/4.0 )
    return static_cast<TFunctionValue>(
        sqrt(0.5 +
          std::tan(this->m_Kappa * (1 + 2 * std::log2(4 * freq_in_hz))) /
          2*std::tan(this->m_Kappa))
        );

  if(freq_in_hz >= 1/4.0  && freq_in_hz <= 1/2.0)
    return static_cast<TFunctionValue>(
        sqrt(0.5 -
          std::tan(this->m_Kappa * (1 + 2 * std::log2(2 * freq_in_hz))) /
          2*std::tan(this->m_Kappa))
        );
  return 0;
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename VowWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
VowWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateForwardLowPassFilter(const FunctionValueType & freq_in_hz) const
{
  FunctionValueType value =
    std::pow(freq_in_hz, this->m_HighPassSubBands) *
    std::pow(2.0, 2*this->m_HighPassSubBands - 1);
  if ( value > 0.25 )
    return this->EvaluateFunction(value);
  return 1;
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename VowWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
VowWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateForwardHighPassFilter(const FunctionValueType & freq_in_hz) const
{
  FunctionValueType value =
    std::pow(freq_in_hz, this->m_HighPassSubBands) *
    std::pow(2.0, this->m_HighPassSubBands - 1);
  if ( value < 0.25 )
    return this->EvaluateFunction(value);
  return 1;
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename VowWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
VowWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateForwardSubBand(const FunctionValueType & freq_in_hz, unsigned int j)
const
{
  if (j == this->m_HighPassSubBands) return this->EvaluateForwardHighPassFilter(freq_in_hz);
  if (j == 0) return this->EvaluateForwardLowPassFilter(freq_in_hz);
  if (j > this->m_HighPassSubBands || j < 0)
    throw itk::ExceptionObject(__FILE__, __LINE__,
            "Invalid SubBand", ITK_LOCATION);
  FunctionValueType value =
    std::pow(freq_in_hz, this->m_HighPassSubBands) *
    std::pow(2.0, 2*this->m_HighPassSubBands - 1 - j);
  return this->EvaluateFunction(value);
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename VowWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
VowWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateInverseLowPassFilter(const FunctionValueType & freq_in_hz) const
{
  return this->EvaluateForwardLowPassFilter(freq_in_hz);
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename VowWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
VowWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateInverseHighPassFilter(const FunctionValueType & freq_in_hz) const
{
  return this->EvaluateForwardHighPassFilter(freq_in_hz);
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename VowWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
VowWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateInverseSubBand(const FunctionValueType & freq_in_hz, unsigned int j) const
{
  return this->EvaluateForwardSubBand(freq_in_hz, j);
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
void
VowWavelet< TFunctionValue, VImageDimension, TInput >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Kappa: " << this->m_Kappa << std::endl;
}
} // end namespace itk

#endif
