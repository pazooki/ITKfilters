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
#ifndef itkHeldWavelet_hxx
#define itkHeldWavelet_hxx

#include <cmath>
#include "itkMath.h"
#include "itkHeldWavelet.h"

namespace itk
{
template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
HeldWavelet< TFunctionValue, VImageDimension, TInput >
::HeldWavelet()
 :m_PolynomialOrder(3)
{}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
HeldWavelet< TFunctionValue, VImageDimension, TInput >
::~HeldWavelet()
{}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename HeldWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
HeldWavelet< TFunctionValue, VImageDimension, TInput >
::Evaluate(const TInput & position) const
{
  TFunctionValue normPosition = 0.0;
  for (unsigned int i = 0; i < VImageDimension ; i++)
    normPosition += position[i]*position[i];

  return this->EvaluateFunction(normPosition);
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename HeldWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
HeldWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateFunction(const FunctionValueType & normPosition) const
{
  if(normPosition <= 0.25 && normPosition > 0.125)
    return static_cast<TFunctionValue>(std::cos(2.0 * itk::Math::pi  *
        this->ComputePolynom(normPosition, this->m_PolynomialOrder)));

  if(normPosition <= 0.5 && normPosition > 0.25)
    return static_cast<TFunctionValue>(std::sin(2.0 * itk::Math::pi  *
        this->ComputePolynom(normPosition / 2.0, this->m_PolynomialOrder)));
  else
    return 0;
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename HeldWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
HeldWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateForwardLowPassFilter(const FunctionValueType & normPosition) const
{
  FunctionValueType value =
    std::pow(normPosition, this->m_HighPassSubBands) *
    std::pow(2.0, 2*this->m_HighPassSubBands - 1);
  if ( value > 0.25 )
    return this->EvaluateFunction(value);
  else
    return 1;
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename HeldWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
HeldWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateForwardHighPassFilter(const FunctionValueType & normPosition) const
{
  FunctionValueType value =
    std::pow(normPosition, this->m_HighPassSubBands) *
    std::pow(2.0, this->m_HighPassSubBands - 1);
  if ( value < 0.25 )
    return this->EvaluateFunction(value);
  else
    return 1;
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename HeldWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
HeldWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateForwardSubBand(const FunctionValueType & normPosition, unsigned int j)
const
{
  if (j == this->m_HighPassSubBands) return this->EvaluateForwardHighPassFilter(normPosition);
  if (j == 0) return this->EvaluateForwardLowPassFilter(normPosition);
  if (j > this->m_HighPassSubBands || j < 0)
    throw itk::ExceptionObject(__FILE__, __LINE__,
            "Invalid SubBand", ITK_LOCATION);
  FunctionValueType value =
    std::pow(normPosition, this->m_HighPassSubBands) *
    std::pow(2.0, 2*this->m_HighPassSubBands - 1 - j);
  return this->EvaluateFunction(value);
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename HeldWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
HeldWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateInverseLowPassFilter(const FunctionValueType & normPosition) const
{
  return this->EvaluateForwardLowPassFilter(normPosition);
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename HeldWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
HeldWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateInverseHighPassFilter(const FunctionValueType & normPosition) const
{
  return this->EvaluateForwardHighPassFilter(normPosition);
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename HeldWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
HeldWavelet< TFunctionValue, VImageDimension, TInput >
::EvaluateInverseSubBand(const FunctionValueType & normPosition, unsigned int j) const
{
  return this->EvaluateForwardSubBand(normPosition, j);
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
typename HeldWavelet< TFunctionValue, VImageDimension, TInput >::FunctionValueType
HeldWavelet< TFunctionValue, VImageDimension, TInput >
::ComputePolynom(const FunctionValueType& x, const unsigned int & order) const
{
  FunctionValueType y = 0.0;
  switch (order) {
  case 0: {
    y = 0.5
      - 2 * x;
    break;
  }
  case 1: {
    y = -1
      + 24  * x
      - 144 * std::pow(x, 2.0)
      + 256 * std::pow(x, 3.0);
    break;
  }
  case 2: {
    y = 8
      - 239   * x
      + 2879  * std::pow(x, 2.0)
      - 16639 * std::pow(x, 3.0)
      + 46079 * std::pow(x, 4.0)
      - 49151 * std::pow(x, 5.0);
    break;
  }
  case 3: {
    y = -52
      + 1120     * x
      - 40320    * std::pow(x, 2.0)
      + 394240   * std::pow(x, 3.0)
      - 2257920  * std::pow(x, 4.0)
      + 7569408  * std::pow(x, 5.0)
      - 13762560 * std::pow(x, 6.0)
      + 10485760 * std::pow(x, 7.0);
    break;
  }
  case 4: {
    FunctionValueType x8 = 8 * x;
    y = 368.0
      - 2520.0  * x8
      + 7560.0  * std::pow(x8, 2.0)
      - 13020.0 * std::pow(x8, 3.0)
      + 14175.0 * std::pow(x8, 4.0)
      - 10111.5 * std::pow(x8, 5.0)
      + 4725.0  * std::pow(x8, 6.0)
      - 1395.0  * std::pow(x8, 7.0)
      + 236.25  * std::pow(x8, 8.0)
      - 17.5    * std::pow(x8, 9.0);
    break;
  }
  case 5: {
    FunctionValueType x8 = 8 * x;
    y = -2656.0
      + 22176.0  * x8
      - 83160.0  * std::pow(x8, 2.0)
      + 184800.0 * std::pow(x8, 3.0)
      - 270270.0 * std::pow(x8, 4.0)
      + 273042.0 * std::pow(x8, 5.0)
      - 194386.5 * std::pow(x8, 6.0)
      + 97515.0  * std::pow(x8, 7.0)
      - 33783.75 * std::pow(x8, 8.0)
      + 7700.0   * std::pow(x8, 9.0)
      - 1039.5   * std::pow(x8, 10.0)
      + 63.0     * std::pow(x8, 11.0);
    break;
  }
  default: {
    throw itk::ExceptionObject(__FILE__, __LINE__,
      "Choose order of polynom less than 6", ITK_LOCATION);
  }
  }
  return y;
}

template< typename TFunctionValue, unsigned int VImageDimension, typename TInput >
void
HeldWavelet< TFunctionValue, VImageDimension, TInput >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "PolynomialOrder: " << this->m_PolynomialOrder << std::endl;
  os << indent << "HighPassSubBands: " << this->m_HighPassSubBands << std::endl;
}
} // end namespace itk

#endif
