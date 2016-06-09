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
#ifndef itkWaveletFunction_h
#define itkWaveletFunction_h

#include "itkSpatialFunction.h"
#include "itkFloatTypes.h"

namespace itk
{
/** \class WaveletFunction
 * Abstract class for WaveletFunction.
 * The interface consists on EvaluateAlgo functions.
 * Where Algo can be:
 *  Function: MotherWavelet response.
 *  ForwardYYY: Forward/Analysis filter response.
 *  InverseYYY: Inverse/Synthesis filter response.
 *  Where YYY can be:
 *   LowPassFilter
 *   HighPassFilter
 *   SubBandFilter
 *
 *  Data member m_HighPassSubBand refers to the numbers of bands in the highpass filter. Default to 1, just one highpass filter, so no band analysis.
 *
 *  If using EvaluateAlgoYYYSubBandFilter(unsigned int k). Implement it like this:
 *  k = 0 return low-pass response.
 *  k = m_HighPassSubBand return high-pass response.
 *  0 < k < m_HighPassSubBand calls the high pass subbands.
 *
 * For an example \sa itkHeldWavelet
 *
 * \ingroup SpatialFunctions
 * \ingroup ITKWavelet
 */
template< typename TFunctionValue = double,
          unsigned int VImageDimension = 3,
          typename TInput = Point< SpacePrecisionType, VImageDimension > >
class WaveletFunction:
  public SpatialFunction< TFunctionValue, VImageDimension, TInput >
{
public:
  /** Standard class typedefs. */
  typedef WaveletFunction                                     Self;
  typedef SpatialFunction< TFunctionValue, VImageDimension, TInput > Superclass;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(WaveletFunction, SpatialFunction);

  /** Input type for the function. */
  typedef typename Superclass::InputType InputType;

  /** Output type for the function. */
  typedef typename Superclass::OutputType FunctionValueType;
  typedef typename Superclass::OutputType OutputType;

  /** Evaluate the function at a given point position. */
  virtual FunctionValueType Evaluate(const TInput & position) const = 0;

  /** Evaluate the function */
  virtual FunctionValueType EvaluateFunction(const TFunctionValue& x) const = 0;

  /**** Forward/Analysis ***/
  /** Evaluate the low filter response. */
  virtual FunctionValueType EvaluateForwardLowPassFilter(const TFunctionValue& x) const = 0;
  /** Evaluate the highfilter response. */
  virtual FunctionValueType EvaluateForwardHighPassFilter(const TFunctionValue& x) const = 0;
  /** Evaluate the sub-band response.
   * j = 0 evaluates LowFilter, j=m_SubBand evaluates HighFilter */
  virtual FunctionValueType EvaluateForwardSubBand( const TFunctionValue& x,
      unsigned int j) const = 0;

  /**** Inverse/Synthesis ***/
  /** Evaluate the low filter response. */
  virtual FunctionValueType EvaluateInverseLowPassFilter(const TFunctionValue& x) const = 0;
  /** Evaluate the highfilter response. */
  virtual FunctionValueType EvaluateInverseHighPassFilter(const TFunctionValue& x) const = 0;
  /** Evaluate the sub-band response.
   * j = 0 evaluates LowFilter, j=m_SubBand evaluates HighFilter */
  virtual FunctionValueType EvaluateInverseSubBand( const TFunctionValue& x,
      unsigned int j) const = 0;

  /** Gets and sets parameters */
  itkSetMacro(HighPassSubBands, unsigned int);
  itkGetConstMacro(HighPassSubBands, unsigned int);

protected:
  WaveletFunction();
  virtual ~WaveletFunction();
  virtual void PrintSelf(std::ostream & os, Indent indent) const;
  /** Number of HighPassSubBands in the high filter decomposition.
   * Default to one HighPass filter (no subbands) */
  unsigned int m_HighPassSubBands;

private:
  WaveletFunction(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWaveletFunction.hxx"
#endif

#endif
