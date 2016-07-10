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
#ifndef itkForwardWaveletFFT_h
#define itkForwardWaveletFFT_h

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include <itkImageConstIterator.h>
#include <complex>
#include <itkFixedArray.h>
#include <itkImageToImageFilter.h>
// #include "itkWaveletFrequencyFilterBankGenerator.h"

namespace itk
{
/** \class ForwardWaveletFFT
 * @brief Wavelet analysis where input is an FFT image.
 * Aim to be Isotropic.
 *
 * \ingroup ITKWavelet
 */
template< class TInputImage, class TOutputImage, class TWaveletFilterBank >
class ForwardWaveletFFT:
  public ImageToImageFilter< TInputImage, TOutputImage>
{
public:
  /** Standard classs typedefs. */
  typedef ForwardWaveletFFT                               Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Inherit types from Superclass. */
  typedef typename Superclass::InputImageType         InputImageType;
  typedef typename Superclass::OutputImageType        OutputImageType;
  typedef typename Superclass::InputImagePointer      InputImagePointer;
  typedef typename Superclass::OutputImagePointer     OutputImagePointer;
  typedef typename Superclass::InputImageConstPointer InputImageConstPointer;

  typedef typename itk::ImageRegionIterator<OutputImageType>     OutputRegionIterator;
  typedef typename itk::ImageRegionConstIterator<InputImageType> InputRegionConstIterator;
  typedef typename OutputImageType::RegionType                   OutputImageRegionType;

  typedef TWaveletFilterBank                                  WaveletFilterBankType;
  typedef typename WaveletFilterBankType::WaveletFunctionType WaveletFunctionType;
  typedef typename WaveletFilterBankType::FunctionValueType   FunctionValueType;

  /** ImageDimension constants */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(ForwardWaveletFFT,
               ImageToImageFilter);
  void SetLevels(unsigned int n);
  itkGetMacro(Levels, unsigned int);
  void SetHighPassSubBands(unsigned int n);
  itkGetMacro(HighPassSubBands, unsigned int);

  /** Compute max number of levels depending on the size of the image.
   * Return J: $ J = min_element(J_0,J_1,...) $ ;
   * where $ 2^J_0 = input_size[0], 2^J_1 = input_size[1] ...  $
   * When the max level is equal to 1 (J=1) implies that there isn't multidimensional analysis.
   * If the size on any dimension is not a power of 2, the max level will be 1.
   * If the sizes are different, but all of them are power of 2, the max level will be the minimum $J_i$.
   */
  unsigned int ComputeMaxNumberOfLevels(typename InputImageType::SizeType& input_size);

protected:
  ForwardWaveletFFT();
  ~ForwardWaveletFFT() {}
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;
  /** Single-threaded version of GenerateData. */
  void GenerateData() ITK_OVERRIDE;

  unsigned int m_Levels;
  unsigned int m_HighPassSubBands;
  unsigned int m_TotalOutputs;
private:
  ForwardWaveletFFT(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
} // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkForwardWaveletFFT.hxx"
#endif

#endif
