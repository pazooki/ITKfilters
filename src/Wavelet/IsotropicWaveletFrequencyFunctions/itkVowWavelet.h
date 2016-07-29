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
#ifndef itkVowWavelet_h
#define itkVowWavelet_h

#include "itkIsotropicWaveletFrequencyFunction.h"

namespace itk
{
/** \class VowWavelet
 * \brief Wavelet based on paper VOW: Variance-Optimal Wavelets for Steerable Pyramid (P.Pad et al 2014).
 *
 * Implement function in frequency space.
 \f{equation}
 \label{VOW}
   h(\omega) =
     \begin{cases}
     \begin{aligned}
       &\sqrt{\frac{1}{2} + \frac{\tan(\kappa(1+2\log_2\frac{2\omega}{\pi})}{2\tan(\kappa)}} , &\omega \in [\frac{\pi}{4} , \frac{\pi}{2} [ \\
       &\sqrt{\frac{1}{2} - \frac{\tan(\kappa(1+2\log_2\frac{\omega}{\pi}))}{2\tan(\kappa)}} , &\omega \in [\frac{\pi}{2} , \pi ] \\
       &0, &\text{otherwise}
     \end{aligned}
     \end{cases}
 \f{equation}
 \f{equation}
   \text{where } \kappa \in [0, \frac{\pi}{2}] \text{ is found to be } 0.75
 \f{equation}
 *
 * Where q(t) is a m grade polynomial (m can be chosen) which elements are
 * calculated so the wavelet has desirable properties.
 * ie, tight frame, Vow Paritition of Unity, etc. (see paper for more info)
 *
 * \ingroup SpatialFunctions
 * \ingroup ITKWavelets
 */
template< typename TFunctionValue = double,
          unsigned int VImageDimension = 3,
          typename TInput = Point< SpacePrecisionType, VImageDimension > >
class VowWavelet:
  public IsotropicWaveletFrequencyFunction< TFunctionValue, VImageDimension, TInput >
{
public:
  /** Standard class typedefs. */
  typedef VowWavelet                                        Self;
  typedef IsotropicWaveletFrequencyFunction< TFunctionValue, VImageDimension, TInput > Superclass;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VowWavelet, SpatialFunction);

  /** Input type for the function. */
  typedef typename Superclass::InputType InputType;

  /** FunctionValue type for the function. */
  typedef typename Superclass::FunctionValueType FunctionValueType;

  /** Type used to store gaussian parameters. */
  typedef FixedArray< double, VImageDimension > ArrayType;

  /** Evaluate the function at a given point position. */
  FunctionValueType Evaluate(const TInput & position) const ITK_OVERRIDE;

  /** Evaluate the function */
  FunctionValueType EvaluateFunction(const FunctionValueType& freq_in_hz) const ITK_OVERRIDE;

  /**** Forward/Analysis ***/
  /** Evaluate the low filter response. */
  FunctionValueType EvaluateForwardLowPassFilter(const FunctionValueType& freq_in_hz) const ITK_OVERRIDE;
  /** Evaluate the highfilter response. */
  FunctionValueType EvaluateForwardHighPassFilter(const FunctionValueType& freq_in_hz) const ITK_OVERRIDE;
  /** Evaluate the sub-band response.
   * j = 0 evaluates LowFilter, j=m_SubBand evaluates HighFilter */
  FunctionValueType EvaluateForwardSubBand( const FunctionValueType& freq_in_hz,
      unsigned int j) const ITK_OVERRIDE;
  /**** Inverse/Synthesis ***/
  /** Evaluate the low filter response. */
  FunctionValueType EvaluateInverseLowPassFilter(const FunctionValueType& freq_in_hz) const ITK_OVERRIDE;
  /** Evaluate the highfilter response. */
  FunctionValueType EvaluateInverseHighPassFilter(const FunctionValueType& freq_in_hz) const ITK_OVERRIDE;
  /** Evaluate the sub-band response.
   * j = 0 evaluates LowFilter, j=m_SubBand evaluates HighFilter */
  FunctionValueType EvaluateInverseSubBand( const FunctionValueType& freq_in_hz,
      unsigned int j) const ITK_OVERRIDE;

  /** Gets and sets parameters */
  itkSetMacro(Kappa, TFunctionValue);
  itkGetConstMacro(Kappa, TFunctionValue);

protected:
  VowWavelet();
  virtual ~VowWavelet();
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:
  VowWavelet(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

  /** kappa value, default is optimal:0.75 */
 TFunctionValue m_Kappa;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVowWavelet.hxx"
#endif

#endif
