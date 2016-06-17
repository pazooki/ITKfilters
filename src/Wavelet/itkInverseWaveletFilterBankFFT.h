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
#ifndef itkInverseWaveletFilterBankFFT_h
#define itkInverseWaveletFilterBankFFT_h

#include <itkImageConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <complex>
#include <itkImageToImageFilter.h>

namespace itk
{

/** \class InverseWaveletFilterBankFFT
 * \brief Low-pass / high-pass wavelet transformation.
 *
 * This implementation performs a low-pass / high-pass wavelet transformation of an image.
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
template <class TInputImage, class TOutputImage, class TWaveletFunction>
class InverseWaveletFilterBankFFT
: public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard typedefs */
  typedef InverseWaveletFilterBankFFT                        Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>                            Pointer;
  typedef itk::SmartPointer<const Self>                      ConstPointer;

  /** Type macro */
  itkNewMacro(Self);

  /** Creation through object factory macro */
  itkTypeMacro(InverseWaveletFilterBankFFT, ImageToImageFilter);

  /** Inherit types from Superclass. */
  typedef typename Superclass::InputImageType         InputImageType;
  typedef typename Superclass::OutputImageType        OutputImageType;
  typedef typename Superclass::InputImagePointer      InputImagePointer;
  typedef typename Superclass::OutputImagePointer     OutputImagePointer;
  typedef typename Superclass::InputImageConstPointer InputImageConstPointer;

  typedef typename itk::ImageRegionIterator<OutputImageType>     OutputRegionIterator;
  typedef typename itk::ImageRegionConstIterator<InputImageType> InputRegionConstIterator;
  typedef typename OutputImageType::RegionType                   OutputImageRegionType;

  typedef TWaveletFunction                                WaveletFunctionType;
  typedef typename WaveletFunctionType::FunctionValueType FunctionValueType;

  /** Dimension */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  /* Getters/Setters */
  /* Members */
  itkGetMacro(ExpandFactor, unsigned int);
  itkGetMacro(HighPassSubBands, unsigned int);
  void SetHighPassSubBands(unsigned int k);
  /* Set Inputs *****/
  void SetInputLowPass(InputImagePointer imgP);
  void SetInputHighPass(InputImagePointer imgP);
  void SetInputSubBand(unsigned int k, InputImagePointer imgP);
  void SetInputs(std::vector<InputImagePointer> inputVector);
protected:
  InverseWaveletFilterBankFFT();
  virtual ~InverseWaveletFilterBankFFT() {}
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

  unsigned int m_ExpandFactor;
  unsigned int m_HighPassSubBands;
private:
  InverseWaveletFilterBankFFT(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self&) ITK_DELETE_FUNCTION;

}; // end of class
} // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkInverseWaveletFilterBankFFT.hxx"
#endif

#endif
