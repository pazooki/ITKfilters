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
#ifndef itkFrequencyShrinkImageFilter_h
#define itkFrequencyShrinkImageFilter_h

#include "itkShrinkImageFilter.h"
#include "itkEnableIf.h"

namespace itk
{

/** \class FrequencyShrinkImageFilter
 * \brief Reduce the size of an image in the FFT domain by an integer factor
 * in each dimension.
 * This filter discard all the high frequency bins depending on the shrink factor
 * | 0        1   ... N/2-1   N/2             : N/2+1 ... N-1 |
 * | 0 (AC)   Low ... High    Highest-Nyquist : High  ... Low |
 * Considering a shrink factor of 2.
 * | 0  1 ... N/4 ... N/2 : N/2+1 ... N-1-N/4 ... N -1 |
 * The outputs will be, no average, just chopping high-frequencies.
 * | 0  1 ... N/4 : N-1-N/4 ... N-1 |
 * If N = Even:
 * | Samples: N/4 + 1 (0 + "New Nyquist") | N/4 - 1 |
 *
 *
 * The output image size in each dimension is given by:
 *
 * outputSize[j] = max( std::floor(inputSize[j]/shrinkFactor[j]), 1 );
 *
 * This filter is implemented so that the starting extent of the first
 * pixel of the output matches that of the input.
 *
 * \image html FrequencyShrinkGrid.png "The change in image geometry from a 5x5 image binned by a factor of 2x2."
 *
 * This code was contributed in the Insight Journal paper:
 * https://hdl.handle.net....
 *
 * \ingroup ITKImageGrid
 */
template <typename TInputImage, typename TOutputImage>
class FrequencyShrinkImageFilter :
  public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef FrequencyShrinkImageFilter                         Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FrequencyShrinkImageFilter, ImageToImageFilter);

  /** Typedef to images */
  typedef TOutputImage                          OutputImageType;
  typedef TInputImage                           InputImageType;
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;

  typedef typename TOutputImage::OffsetType  OutputOffsetType;
  typedef typename TOutputImage::IndexType   OutputIndexType;
  typedef typename TInputImage::IndexType    InputIndexType;

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension );
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension );

  typedef FixedArray< unsigned int, ImageDimension > ShrinkFactorsType;

  /** Set the shrink factors. Values are clamped to
   * a minimum value of 1. Default is 1 for all dimensions. */
  itkSetMacro(ShrinkFactors, ShrinkFactorsType);
  void SetShrinkFactors(unsigned int factor);
  void SetShrinkFactor(unsigned int i, unsigned int factor);

  /** Get the shrink factors. */
  itkGetConstReferenceMacro(ShrinkFactors, ShrinkFactorsType);

  virtual void GenerateOutputInformation() ITK_OVERRIDE;

  /** FrequencyShrinkImageFilter needs a larger input requested region than the output
   * requested region.  As such, FrequencyShrinkImageFilter needs to provide an
   * implementation for GenerateInputRequestedRegion() in order to inform the
   * pipeline execution model.
   * \sa ProcessObject::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;


#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToOutputCheck,
    (Concept::Convertible<typename TInputImage::PixelType, typename TOutputImage::PixelType>));
  itkConceptMacro(SameDimensionCheck,
    (Concept::SameDimension<ImageDimension, OutputImageDimension>));
  /** End concept checking */
#endif

protected:
  FrequencyShrinkImageFilter();
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  void GenerateData() ITK_OVERRIDE;

private:
  FrequencyShrinkImageFilter(const Self&) ITK_DELETE_FUNCTION;
  void operator=(const Self&) ITK_DELETE_FUNCTION;

  ShrinkFactorsType m_ShrinkFactors;

  /** Round different pixel types. */
  template< class TOutputType, class TInputType >
  typename EnableIfC<std::numeric_limits<TOutputType>::is_integer,  TOutputType>::Type
  RoundIfInteger( TInputType input )
    {
      return Math::Round< TOutputType >( input );
    }

  // For Non-fundamental types numeric_limits is not specialized, and
  // is_integer defaults to false.
  template< class TOutputType, class TInputType >
  typename DisableIfC<std::numeric_limits<TOutputType>::is_integer,  TOutputType>::Type
  RoundIfInteger( const TInputType & input, ...)
    {
      return static_cast<TOutputType>(input);
    }
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFrequencyShrinkImageFilter.hxx"
#endif

#endif
