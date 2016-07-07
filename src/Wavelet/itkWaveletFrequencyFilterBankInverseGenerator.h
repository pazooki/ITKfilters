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
#ifndef itkWaveletFrequencyFilterBankInverseGenerator_h
#define itkWaveletFrequencyFilterBankInverseGenerator_h

#include <itkImageConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <complex>
#include <itkGenerateImageSource.h>

namespace itk
{

/** \class WaveletFrequencyFilterBankInverseGenerator
 * \brief Generate filter bank of inverse wavelet images.
 *
 * This implementation generates a low-pass / high-pass wavelet filter .
 *
 * The output are two images, one low-pass and other high-pass, downsampled by DownSampleImageFactor.
 *
 * The wavelet operation is defined by an user chosen SpatialFunction
 * The Spatial Function can be in the frequency or spatial domain.
 *
 * \sa WaveletFT
 * \sa SpatialFunction
 *
 * \ingroup ITKWavelet
 */
template < class TOutputImage, class TWaveletFunction>
class WaveletFrequencyFilterBankInverseGenerator
: public itk::GenerateImageSource< TOutputImage >
{
public:
  /** Standard typedefs */
  typedef WaveletFrequencyFilterBankInverseGenerator            Self;
  typedef itk::GenerateImageSource<TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>                Pointer;
  typedef itk::SmartPointer<const Self>          ConstPointer;

  /** Type macro */
  itkNewMacro(Self);

  /** Creation through object factory macro */
  itkTypeMacro(WaveletFrequencyFilterBankInverseGenerator, GenerateImageSource);

  /** Inherit types from Superclass. */
  typedef typename Superclass::OutputImageType        OutputImageType;
  typedef typename Superclass::OutputImagePointer     OutputImagePointer;

  typedef typename itk::ImageRegionIterator<OutputImageType>     OutputRegionIterator;
  typedef typename OutputImageType::RegionType                   OutputImageRegionType;

  typedef TWaveletFunction                                WaveletFunctionType;
  typedef typename WaveletFunctionType::FunctionValueType FunctionValueType;

  /** Dimension */
  itkStaticConstMacro(ImageDimension, unsigned int, TOutputImage::ImageDimension);

  /* Getters/Setters */
  /* Members */
  itkGetMacro(HighPassSubBands, unsigned int);
  void SetHighPassSubBands(unsigned int k);
  itkGetMacro(ScaleFactor, unsigned int);
  itkSetMacro(ScaleFactor, unsigned int);
  /** Get Outputs *****/
  OutputImagePointer GetOutputLowPass();
  OutputImagePointer GetOutputHighPass();
  OutputImagePointer GetOutputSubBand(unsigned int k);
  std::vector<OutputImagePointer> GetOutputs();
protected:
  WaveletFrequencyFilterBankInverseGenerator();
  virtual ~WaveletFrequencyFilterBankInverseGenerator() {}
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  /************ GenerateData *************/

  /** BeforeThreadedGenerateData.
   * It allocates also internal images
   */
  // virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;

  /** AfterThreadedGenerateData.
   * It enforce memory destruction of internal images
   */
  // virtual void AfterThreadedGenerateData() ITK_OVERRIDE;

  /** Generate data redefinition */
  virtual void GenerateData() ITK_OVERRIDE;
  // virtual void ThreadedGenerateData(
  //     const OutputImageRegionType& outputRegionForThread,
  //     itk::ThreadIdType threadId) ITK_OVERRIDE;
  /************ Data Members *************/

  /** Number of M-Bands decomposition of the high pass filters */
  unsigned int m_HighPassSubBands;
  /** Scale Factor, defaults to 1, used optionally to dilate the filter bank in
   * a multiresolution (multi-level) analysis. */
  unsigned int m_ScaleFactor;
private:
  WaveletFrequencyFilterBankInverseGenerator(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self&) ITK_DELETE_FUNCTION;

}; // end of class
} // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWaveletFrequencyFilterBankInverseGenerator.hxx"
#endif

#endif
