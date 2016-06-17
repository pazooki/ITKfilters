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

#include <itkForwardWaveletFilterBankFFT.h>
#include "itkImageRegionConstIterator.h"
#include <itkImageConstIterator.h>
#include <complex>
#include <itkFixedArray.h>

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
  typedef itk::ForwardWaveletFilterBankFFT<OutputImageType, OutputImageType, WaveletFunctionType> OutputWaveletFilterBankType;

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

protected:
  ForwardWaveletFFT();
  ~ForwardWaveletFFT() {}
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;
  /** Single-threaded version of GenerateData. */
  void GenerateData() ITK_OVERRIDE;

  unsigned int m_Levels;
  unsigned int m_HighPassSubBands;
  unsigned int m_TotalOutputs;
  typename WaveletFilterBankType::Pointer m_FilterBank;
private:
  ForwardWaveletFFT(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
} // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkForwardWaveletFFT.hxx"
#endif

#endif
