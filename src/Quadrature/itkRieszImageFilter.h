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
#ifndef itkRieszImageFilter_h
#define itkRieszImageFilter_h

#include "itkImageToImageFilter.h"
#include <complex>

namespace itk
{
/** \class RieszImageFilter
 * @brief Analytical filter, analogous to Hilber transform in nD.
 *
 * \ingroup ITKBasicFilters
 */
template< typename TInputImage, typename TOutputImage =
  Image<std::complex<typename TInputImage::PixelType>,
        TInputImage::ImageDimension > >
class RieszImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef RieszImageFilter                                Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                            InputImageType;
  typedef typename InputImageType::Pointer       InputImagePointer;
  typedef typename InputImageType::ConstPointer  InputImageConstPointer;
  typedef typename InputImageType::RegionType    InputImageRegionType;
  typedef typename InputImageType::PixelType     InputImagePixelType;
  // Intermediate calculations:
  typedef Image<std::complex<typename InputImagePixelType>,
          TInputImage::ImageDimension>           FFTImageType;
  typedef typename FFTImageType::Pointer      FFTImagePointer;
  typedef typename FFTImageType::ConstPointer FFTImageConstPointer;
  typedef typename FFTImageType::RegionType   FFTImageRegionType;
  typedef typename FFTImageType::PixelType    FFTImagePixelType;

  /** ImageDimension constants */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(RieszImageFilter,
               ImageToImageFilter);
protected:
  RieszImageFilter();
  ~RieszImageFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;
  /** Single-threaded version of GenerateData. */
  void GenerateData() ITK_OVERRIDE;
private:
  RieszImageFilter(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
} // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRieszImageFilter.hxx"
#endif

#endif
