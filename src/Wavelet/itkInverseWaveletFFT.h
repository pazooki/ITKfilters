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
#ifndef itkInverseWaveletFT_h
#define itkInverseWaveletFT_h

#include <itkImageToImageListFilter.h>
#include <itkVectorImage.h>
#include "itkStatisticsImageFilter.h"
#include <itkImageConstIterator.h>
#include <complex>
#include <itkFixedArray.h>
#include <itkSymmetricSecondRankTensor.h>
#include "itkRealToHalfHermitianForwardFFTImageFilter.h"
#include "itkHalfHermitianToRealInverseFFTImageFilter.h"
// #ifdef WITH_C++11
#include <functional>
// #endif

namespace itk
{
/** \class InverseWaveletFT
 * @brief Analytical filter, analogous to Hilber transform in nD.
 *
 * \ingroup ITKBasicFilters
 */
template< typename TInputImage >
class InverseWaveletFT:
  public ImageToImageListFilter< TInputImage, TInputImage>
{
public:
  /** Standard class typedefs. */
  typedef InverseWaveletFT                                Self;
  typedef ImageToImageFilter< TInputImage, TInputImage >  Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  typedef typename InputImageType::PixelType      InputImagePixelType;
  typedef typename InputImageType::SpacingType    SpacingType;
  typedef typename InputImageRegionType::SizeType SizeType;
  // Intermediate calculations:
  typedef Image<std::complex<InputImagePixelType>,
          TInputImage::ImageDimension>            ComplexImageType;
  typedef typename ComplexImageType::Pointer      ComplexImagePointer;
  typedef typename ComplexImageType::ConstPointer ComplexImageConstPointer;
  typedef typename ComplexImageType::RegionType   ComplexImageRegionType;
  typedef typename ComplexImageType::PixelType    ComplexImagePixelType;

  typedef itk::StatisticsImageFilter<InputImageType> StatisticsImageFilterType;
  typedef typename StatisticsImageFilterType::RealType StatisticsRealType;

  typedef double RealType;
  typedef FixedArray< RealType, TInputImage::ImageDimension > DirectionType;
  typedef VectorImage< InputImagePixelType, TInputImage::ImageDimension >
    RieszComponentsImageType;
  typedef VectorImage< ComplexImagePixelType, TInputImage::ImageDimension >
    ComplexRieszComponentsImageType;
  typedef RealToHalfHermitianForwardFFTImageFilter < InputImageType>
    FFTFilterType;
  typedef HalfHermitianToRealInverseFFTImageFilter < ComplexImageType,
          InputImageType> InverseFFTFilterType;

  /** ImageDimension constants */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(InverseWaveletFT,
               ImageToImageFilter);

protected:
  InverseWaveletFT();
  ~InverseWaveletFT() {}
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;
  /** Single-threaded version of GenerateData. */
  void GenerateData() ITK_OVERRIDE;
  using Superclass::MakeOutput;
  typedef ProcessObject::DataObjectPointerArraySizeType
    DataObjectPointerArraySizeType; /**  Create the Outputs */
  virtual DataObject::Pointer MakeOutput(
      DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

private:
  InverseWaveletFT(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
} // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkInverseWaveletFT.hxx"
#endif

#endif
