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

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImageConstIterator.h>
#include <complex>
#include <itkFixedArray.h>
#include <itkSymmetricSecondRankTensor.h>

namespace itk
{
/** \class RieszImageFilter
 * @brief Analytical filter, analogous to Hilber transform in nD.
 *
 * \ingroup ITKBasicFilters
 */
template< typename TInputImage >
class RieszImageFilter:
  public ImageToImageFilter< TInputImage, TInputImage>
{
public:
  /** Standard class typedefs. */
  typedef RieszImageFilter                                Self;
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

  typedef double RealType;
  typedef FixedArray< RealType, TInputImage::ImageDimension > DirectionType;
  typedef VectorImage< InputImagePixelType, TInputImage::ImageDimension > RieszComponentsImageType;
  typedef VectorImage< ComplexImagePixelType, TInputImage::ImageDimension > ComplexRieszComponentsImageType;

  // typedef Image< FixedArray<InputImagePixelType, TInputImage::ImageDimension >, TInputImage::ImageDimension > RieszComponentsArrayImageType;
  /** ImageDimension constants */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(RieszImageFilter,
               ImageToImageFilter);

  itkSetMacro( SigmaGaussianDerivative, RealType );
  itkGetConstMacro( SigmaGaussianDerivative, RealType );


  TInputImage*              GetOutputReal();
  RieszComponentsImageType* GetOutputRieszComponents();
  TInputImage*              GetOutputRieszNorm();

  typename InputImageType::Pointer ComputeRieszProjection(const DirectionType & direction, const RieszComponentsImageType* rieszComponents ) const;
  RealType ComputeLocalRieszProjection(const DirectionType & direction, const ImageConstIterator<RieszComponentsImageType> & rieszIt ) const;
  typename InputImageType::Pointer ComputeLocalPhaseInDirection(const DirectionType & unitary_direction, const InputImageType* rieszReal, const RieszComponentsImageType* rieszComponents) const;

  typedef itk::Image<itk::Matrix<RealType,ImageDimension, ImageDimension>, ImageDimension >
    EigenVectorsImageType;
  typedef itk::Image<itk::FixedArray<RealType,ImageDimension>, ImageDimension > EigenValuesImageType;

  std::pair<typename EigenVectorsImageType::Pointer, typename EigenValuesImageType::Pointer>
    ComputeEigenAnalysisMaximizingRieszComponents(const unsigned int & gaussian_window_radius, const float & gaussian_window_sigma, const RieszComponentsImageType* rieszComponents) const;

  typename RieszComponentsImageType::Pointer
    ComputeRieszComponentsWithMaximumResponse(const typename EigenVectorsImageType::Pointer eigenVectors, const RieszComponentsImageType* rieszComponents) const;

  typename InputImageType::Pointer ComputeRealComponent(const ComplexImageType* fftForward) const;

  typename InputImageType::Pointer ComputeLocalAmplitude(
        const InputImageType* real_part,
        const InputImageType* riesz_norm_part) const;

  typename InputImageType::Pointer ComputeRieszWeightedNorm(
      const RieszComponentsImageType* rieszComponents,
      const DirectionType & weights) const;

  typename InputImageType::Pointer ComputeRieszNorm(
      const RieszComponentsImageType* rieszComponents) const;

  typename InputImageType::Pointer ComputeRieszWeightedNormByEigenValues(
      const RieszComponentsImageType* rieszComponents,
      const typename EigenValuesImageType::Pointer eigenValues) const;
protected:
  RieszImageFilter();
  ~RieszImageFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;
  /** Single-threaded version of GenerateData. */
  void GenerateData() ITK_OVERRIDE;
  using Superclass::MakeOutput;
  typedef ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType; /**  Create the Outputs */
  virtual DataObject::Pointer MakeOutput(DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

private:
  RieszImageFilter(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
  InputImagePixelType m_SigmaGaussianDerivative;
  SizeType       m_inputSizeSquare;    //size of input/output image
  SpacingType    m_inputSpacingSquare; //spacing of input image.

  typename InputImageType::Pointer ComputeRieszComponent(const ComplexImageType* fftForward, const unsigned int & NComponent) const;
  typename RieszComponentsImageType::Pointer ComputeRieszComponents(const ComplexImageType* fftForward) const;
};
} // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRieszImageFilter.hxx"
#endif

#endif
