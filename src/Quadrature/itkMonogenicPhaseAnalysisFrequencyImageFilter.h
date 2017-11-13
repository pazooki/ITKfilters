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
#ifndef itkMonogenicPhaseAnalysisFrequencyImageFilter_h
#define itkMonogenicPhaseAnalysisFrequencyImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImageRegionConstIterator.h>
#include <itkImage.h>
#include <itkFixedArray.h>
#include <itkSymmetricSecondRankTensor.h>

#include <itkFrequencyImageRegionIteratorWithIndex.h>
namespace itk
{
/** \class MonogenicPhaseAnalysisFrequencyImageFilter
 * Analytical filter, analogous to Hilbert transform in nD.
 * Require input is a complex image.
 *
 * \ingroup IsotropicWavelets
 */
template< typename TInputImage,
          typename TFrequencyImageRegionConstIterator =
            FrequencyImageRegionIteratorWithIndex< TInputImage> > >
class MonogenicPhaseAnalysisFrequencyImageFilter:
  public ImageToImageFilter<TInputImage, TInputImage>
    {
public:
  /** Standard class typedefs. */
  typedef MonogenicPhaseAnalysisFrequencyImageFilter                                Self;
  typedef ImageToImageFilter< TInputImage, TInputImage >  Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** ImageDimension constants */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(MonogenicPhaseAnalysisFrequencyImageFilter,
               ImageToImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
 /// This ensure that InputPixelType is complex<float||double>
  itkConceptMacro( InputPixelTypeIsComplexAndFloatCheck,
                   ( Concept::IsFloatingPoint< typename TInputImage::PixelType::value_type > ) );
#endif

  /** Some convenient typedefs. */
  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef TFrequencyRegionIterator             OutputFrequencyRegionIterator;
  typedef OutputImageType                      RieszComponentsImageType;

  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  typedef typename InputImageType::PixelType      InputImagePixelType;
  typedef typename InputImageType::SpacingType    SpacingType;
  typedef typename InputImageRegionType::SizeType SizeType;

  typedef SpacingType                               DirectionType;
  typedef typename InputImageType::SpacingValueType FloatType;

  typedef typename OutputImageType::Pointer                       OutputImagePointer;
  typedef typename OutputImageType::ConstPointer                  OutputImageConstPointer;
  typedef typename itk::ImageRegionConstIterator<OutputImageType> OutputImageRegionConstIterator;
  typename InputImageType::Pointer ComputeRieszProjection(
      const DirectionType & direction) const;
      // const RieszComponentsImageType* rieszComponents ) const;

  FloatType ComputeLocalRieszProjection(
      const DirectionType & direction) const;
      // const ImageConstIterator<RieszComponentsImageType> & rieszIt ) const;

  typename InputImageType::Pointer ComputeLocalPhaseInDirection(
      const DirectionType & unitary_direction) const ;
      // const InputImageType* rieszReal,
      // const RieszComponentsImageType* rieszComponents) const;

  typedef itk::Image<itk::Matrix<FloatType,ImageDimension, ImageDimension>,
          ImageDimension > EigenVectorsImageType;
  typedef itk::Image<itk::FixedArray<FloatType,ImageDimension>,
          ImageDimension > EigenValuesImageType;

  std::pair<
        typename EigenVectorsImageType::Pointer,
        typename EigenValuesImageType::Pointer >
    ComputeEigenAnalysisMaximizingRieszComponents(
        const unsigned int & gaussian_window_radius,
        const float & gaussian_window_sigma) const;
        // const RieszComponentsImageType* rieszComponents) const;

  typename RieszComponentsImageType::Pointer
    ComputeRieszComponentsWithMaximumResponse(
        const typename EigenVectorsImageType::Pointer eigenVectors,
        const RieszComponentsImageType* rieszComponents) const;

  // THis shit is input image multiplied by exp(sigma^2 * w^2).
  // I think this is the Mono_0 (but multiplied by a wavelet-like function)
  // This should be eqv to this->GetInput()
  // typename InputImageType::Pointer ComputeRealComponent(
  //     const ComplexImageType* fftForward) const;

  typename InputImageType::Pointer ComputeRieszNorm(
      const RieszComponentsImageType* rieszComponents) const;

  // This does sqrt(M_0^2 + M_R^2)
  // where M_R is the norm of the Riesz Vector.
  typename InputImageType::Pointer ComputeLocalAmplitude(
        const InputImageType* real_part,
        const InputImageType* riesz_norm_part) const;

  typename InputImageType::Pointer ComputeRieszWeightedNorm(
      const RieszComponentsImageType* rieszComponents,
      const DirectionType & weights) const;

  typename InputImageType::Pointer ComputeRieszWeightedNormByEigenValues(
      const RieszComponentsImageType* rieszComponents,
      const typename EigenValuesImageType::Pointer eigenValues) const;
protected:
  MonogenicPhaseAnalysisFrequencyImageFilter();
  ~MonogenicPhaseAnalysisFrequencyImageFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  virtual void GenerateOutputInformation() ITK_OVERRIDE;
  void GenerateData() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MonogenicPhaseAnalysisFrequencyImageFilter);
};
} // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMonogenicPhaseAnalysisFrequencyImageFilter.hxx"
#endif

#endif
