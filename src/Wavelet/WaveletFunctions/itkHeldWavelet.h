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
#ifndef itkHeldWavelet_h
#define itkHeldWavelet_h

#include "itkWaveletFunction.h"

namespace itk
{
/** \class HeldWavelet
 * \brief Wavelet based on paper Steerable Wavelet Frames Based on the Held
 * Transform (Held et al 2010).
 *
 * Implement function in frequency space.
 *
 * h(w) = cos(2*pi*q(|w|)) for w in (1/8, 1/4]
 * h(w) = sin(2*pi*q(|w/2|)) for w in (1/4, 1/2]
 * h(w) = 0 elsewhere.
 *
 * Where q(t) is a m grade polynomial (m can be chosen) which elements are
 * calculated so the wavelet has desirable properties.
 * ie, tight frame, Held Paritition of Unity, etc. (see paper for more info)
 *
 * \ingroup SpatialFunctions
 * \ingroup ITKCommon
 */
template< typename TFunctionValue = double,
          unsigned int VImageDimension = 3,
          typename TInput = Point< SpacePrecisionType, VImageDimension > >
class HeldWavelet:
  public WaveletFunction< TFunctionValue, VImageDimension, TInput >
{
public:
  /** Standard class typedefs. */
  typedef HeldWavelet                                        Self;
  typedef WaveletFunction< TFunctionValue, VImageDimension, TInput > Superclass;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(HeldWavelet, SpatialFunction);

  /** Input type for the function. */
  typedef typename Superclass::InputType InputType;

  /** FunctionValue type for the function. */
  typedef typename Superclass::FunctionValueType FunctionValueType;

  /** Type used to store gaussian parameters. */
  typedef FixedArray< double, VImageDimension > ArrayType;

  /** Evaluate the function at a given point position. */
  FunctionValueType Evaluate(const TInput & position) const ITK_OVERRIDE;

  /** Evaluate the function */
  FunctionValueType EvaluateFunction(const FunctionValueType& normPosition) const ITK_OVERRIDE;

  /**** Forward/Analysis ***/
  /** Evaluate the low filter response. */
  FunctionValueType EvaluateForwardLowPassFilter(const FunctionValueType& normPosition) const ITK_OVERRIDE;
  /** Evaluate the highfilter response. */
  FunctionValueType EvaluateForwardHighPassFilter(const FunctionValueType& normPosition) const ITK_OVERRIDE;
  /** Evaluate the sub-band response.
   * j = 0 evaluates LowFilter, j=m_SubBand evaluates HighFilter */
  FunctionValueType EvaluateForwardSubBand( const FunctionValueType& normPosition,
      unsigned int j) const ITK_OVERRIDE;
  /**** Inverse/Synthesis ***/
  /** Evaluate the low filter response. */
  FunctionValueType EvaluateInverseLowPassFilter(const FunctionValueType& normPosition) const ITK_OVERRIDE;
  /** Evaluate the highfilter response. */
  FunctionValueType EvaluateInverseHighPassFilter(const FunctionValueType& normPosition) const ITK_OVERRIDE;
  /** Evaluate the sub-band response.
   * j = 0 evaluates LowFilter, j=m_SubBand evaluates HighFilter */
  FunctionValueType EvaluateInverseSubBand( const FunctionValueType& normPosition,
      unsigned int j) const ITK_OVERRIDE;

  /** Gets and sets parameters */
  itkSetMacro(PolynomialOrder, unsigned int);
  itkGetConstMacro(PolynomialOrder, unsigned int);

  FunctionValueType ComputePolynom(
      const FunctionValueType & normPosition,
      const unsigned int & order) const;

protected:
  HeldWavelet();
  virtual ~HeldWavelet();
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:
  HeldWavelet(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

  /** The order of the polynom. */
  unsigned int m_PolynomialOrder;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHeldWavelet.hxx"
#endif

#endif
