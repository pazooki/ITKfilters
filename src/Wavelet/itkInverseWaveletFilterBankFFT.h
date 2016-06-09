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
#ifndef itkInverseWaveletFTFilterBank_h
#define itkInverseWaveletFTFilterBank_h

#include <itkImageToImageListFilter.h>
#include <itkVectorImage.h>
#include "itkStatisticsImageFilter.h"
#include <itkImageConstIterator.h>
#include <complex>
#include <itkFixedArray.h>
#include <itkSymmetricSecondRankTensor.h>

namespace itk
{

/** \class InverseWaveletFTFilterBank
 * \brief Low-pass / high-pass wavelet transformation.
 *
 * This implementation performs an inverse Wavelet transformation from a low-pass / high-pass images.
 *
 * The output is the synthesis of the input low-pass/high-pass images.
 *
 * The wavelet operation is defined by an user chosen SpatialFunction
 * The Spatial Function can be in the frequency or spatial domain.
 *
 * \sa WaveletFTFilterBank
 *
 * \ingroup ITKWavelet
 */
template <class TInputImage, class TOutputImage, class TWaveletFunction>
class InverseWaveletFTFilterBank
: public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard typedefs */
  typedef InverseWaveletFTFilterBank                                Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>                            Pointer;
  typedef itk::SmartPointer<const Self>                      ConstPointer;

  /** Type macro */
  itkNewMacro(Self);

  /** Creation through object factory macro */
  itkTypeMacro(InverseWaveletFTFilterBank, ImageToImageFilter);

  /** Template parameters typedefs */
  typedef TInputImage                          InputImageType;
  typedef typename InputImageType::Pointer     InputImagePointerType;
  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename InputImageType::SizeType    InputSizeType;
  typedef typename InputImageType::IndexType   InputIndexType;
  typedef typename InputImageType::PixelType   InputPixelType;

  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointerType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;
  typedef typename OutputImageType::SizeType   OutputSizeType;
  typedef typename OutputImageType::IndexType  OutputIndexType;
  typedef typename OutputImageType::PixelType  OutputPixelType;

  typedef TWaveletFunction                     WaveletFunctionType;

  /** Dimension */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  itkGetMacro(UpSampleImageFactor, unsigned int);
  itkSetMacro(UpSampleImageFactor, unsigned int);

protected:
  InverseWaveletFTFilterBank();
  virtual ~InverseWaveletFTFilterBank() {}
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  /** InverseWaveletFTFilterBank produces images which are of
   * different resolution and different pixel spacing than its input image.
   * As such, InverseWaveletFTFilterBank needs to provide an
   * implementation for GenerateOutputInformation() in order to inform the
   * pipeline execution model.  The original documentation of this method is
   * below.
   * \sa ProcessObject::GenerateOutputInformaton()
   */
  virtual void GenerateOutputInformation() ITK_OVERRIDE;

  /** Given one output whose requested region has been set, this method sets
   * the requested region for the remaining output images.  The original
   * documentation of this method is below.
   * \sa ProcessObject::GenerateOutputRequestedRegion()
   */
  virtual void GenerateOutputRequestedRegion(DataObject *output) ITK_OVERRIDE;

  /** InverseWaveletFTFilterBank requires a larger input requested
   * region than the output requested regions to accommodate the shrinkage and
   * smoothing operations. As such, InverseWaveletFTFilterBank needs
   * to provide an implementation for GenerateInputRequestedRegion().  The
   * original documentation of this method is below.
   * \sa ProcessObject::GenerateInputRequestedRegion()
   */
  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

  /** BeforeThreadedGenerateData.
   * It allocates also internal images
   */
  virtual void BeforeThreadedGenerateData();

  /** Generate data redefinition */
  virtual void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, itk::ThreadIdType threadId);

  /** AfterThreadedGenerateData.
   * It enforce memory destruction of internal images
   */
  virtual void AfterThreadedGenerateData();

  unsigned int m_UpSampleImageFactor;

private:
  InverseWaveletFTFilterBank(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self&) ITK_DELETE_FUNCTION;

}; // end of class
} // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkInverseWaveletFTFilterBank.hxx"
#endif

#endif
