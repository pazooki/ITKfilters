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
#ifndef itkInverseWaveletFTFilterBank_hxx
#define itkInverseWaveletFTFilterBank_hxx
#include "itkInverseWaveletFTFilterBank.h"

namespace itk
{
template< typename TInputImage, typename TOutputImage, typename TWaveletFunction>
InverseWaveletFTFilterBank< TInputImage, TOutputImage, TWaveletFunction>
::InverseWaveletFTFilterBank()
{
  this->SetNumberOfRequiredInputs(2);
  this->SetNumberOfRequiredOutputs(1);

  this->SetNthOutput( 0, this->MakeOutput(0) );

  this->m_UpSampleFilterFactorImageFactor = 2;

}

template< typename TInputImage, typename TOutputImage, typename TWaveletFunction>
void
InverseWaveletFTFilterBank< TInputImage, TOutputImage, TWaveletFunction>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "UpSampleImageFactor : " << m_UpSampleImageFactor
     << std::endl;
}

/*
 * GenerateOutputInformation
 */
template< typename TInputImage, typename TOutputImage, typename TWaveletFunction >
void
InverseWaveletFTFilterBank< TInputImage, TOutputImage, TWaveletFunction >
::GenerateOutputInformation()
{
  // call the superclass's implementation of this method
  Superclass::GenerateOutputInformation();

  // get pointers to the input and output
  InputImageConstPointer inputPtr = this->GetInput();

  if ( !inputPtr  )
    {
    itkExceptionMacro(<< "Input has not been set");
    }

  const typename InputImageType::PointType &
  inputOrigin = inputPtr->GetOrigin();
  const typename InputImageType::SpacingType &
  inputSpacing = inputPtr->GetSpacing();
  const typename InputImageType::DirectionType &
  inputDirection = inputPtr->GetDirection();
  const typename InputImageType::SizeType & inputSize =
    inputPtr->GetLargestPossibleRegion().GetSize();
  const typename InputImageType::IndexType & inputStartIndex =
    inputPtr->GetLargestPossibleRegion().GetIndex();

  typedef typename OutputImageType::SizeType  SizeType;
  typedef typename OutputImageType::IndexType IndexType;

  OutputImagePointer outputPtr;
  typename OutputImageType::PointType outputOrigin;
  typename OutputImageType::SpacingType outputSpacing;
  SizeType  outputSize;
  IndexType outputStartIndex;

  // we need to compute the output spacing, the output image size,
  // and the output image start index
  for ( unsigned int ilevel = 0; ilevel < m_NumberOfLevels; ilevel++ )
    {
    outputPtr = this->GetOutput(ilevel);
    if ( !outputPtr ) { continue; }

    for ( unsigned int idim = 0; idim < OutputImageType::ImageDimension; idim++ )
      {
      const double shrinkFactor = static_cast< double >( m_Schedule[ilevel][idim] );
      outputSpacing[idim] = inputSpacing[idim] * shrinkFactor;

      outputSize[idim] = static_cast< SizeValueType >(
        std::floor(static_cast< double >( inputSize[idim] ) / shrinkFactor) );
      if ( outputSize[idim] < 1 ) { outputSize[idim] = 1; }

      outputStartIndex[idim] = static_cast< IndexValueType >(
        std::ceil(static_cast< double >( inputStartIndex[idim] ) / shrinkFactor) );
      }
    //Now compute the new shifted origin for the updated levels;
    const typename OutputImageType::PointType::VectorType outputOriginOffset =
      ( inputDirection * ( outputSpacing - inputSpacing ) ) * 0.5;
    for ( unsigned int idim = 0; idim < OutputImageType::ImageDimension; idim++ )
      {
      outputOrigin[idim] = inputOrigin[idim] + outputOriginOffset[idim];
      }

    typename OutputImageType::RegionType outputLargestPossibleRegion;
    outputLargestPossibleRegion.SetSize(outputSize);
    outputLargestPossibleRegion.SetIndex(outputStartIndex);

    outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);
    outputPtr->SetOrigin (outputOrigin);
    outputPtr->SetSpacing(outputSpacing);
    outputPtr->SetDirection(inputDirection);  //Output Direction should be same
                                              // as input.
    }
}

/*
 * GenerateOutputRequestedRegion
 */
template< typename TInputImage, typename TOutputImage, typename TWaveletFunction >
void
InverseWaveletFTFilterBank< TInputImage, TOutputImage, TWaveletFunction >
::GenerateOutputRequestedRegion(DataObject *refOutput)
{
  // call the superclass's implementation of this method
  Superclass::GenerateOutputRequestedRegion(refOutput);

  // find the index for this output
  unsigned int refLevel = static_cast<unsigned int>( refOutput->GetSourceOutputIndex() );

  // compute baseIndex and baseSize
  typedef typename OutputImageType::SizeType   SizeType;
  typedef typename OutputImageType::IndexType  IndexType;
  typedef typename OutputImageType::RegionType RegionType;

  TOutputImage *ptr = itkDynamicCastInDebugMode< TOutputImage * >( refOutput );
  if ( !ptr )
    {
    itkExceptionMacro(<< "Could not cast refOutput to TOutputImage*.");
    }

  unsigned int ilevel, idim;

  if ( ptr->GetRequestedRegion() == ptr->GetLargestPossibleRegion() )
    {
    // set the requested regions for the other outputs to their
    // requested region

    for ( ilevel = 0; ilevel < m_NumberOfLevels; ilevel++ )
      {
      if ( ilevel == refLevel ) { continue; }
      if ( !this->GetOutput(ilevel) ) { continue; }
      this->GetOutput(ilevel)->SetRequestedRegionToLargestPossibleRegion();
      }
    }
  else
    {
    // compute requested regions for the other outputs based on
    // the requested region of the reference output
    IndexType  outputIndex;
    SizeType   outputSize;
    RegionType outputRegion;
    IndexType  baseIndex = ptr->GetRequestedRegion().GetIndex();
    SizeType   baseSize  = ptr->GetRequestedRegion().GetSize();

    for ( idim = 0; idim < TOutputImage::ImageDimension; idim++ )
      {
      unsigned int factor = m_Schedule[refLevel][idim];
      baseIndex[idim] *= static_cast< IndexValueType >( factor );
      baseSize[idim] *= static_cast< SizeValueType >( factor );
      }

    for ( ilevel = 0; ilevel < m_NumberOfLevels; ilevel++ )
      {
      if ( ilevel == refLevel ) { continue; }
      if ( !this->GetOutput(ilevel) ) { continue; }

      for ( idim = 0; idim < TOutputImage::ImageDimension; idim++ )
        {
        double factor = static_cast< double >( m_Schedule[ilevel][idim] );

        outputSize[idim] = static_cast< SizeValueType >(
          std::floor(static_cast< double >( baseSize[idim] ) / factor) );
        if ( outputSize[idim] < 1 ) { outputSize[idim] = 1; }

        outputIndex[idim] = static_cast< IndexValueType >(
          std::ceil(static_cast< double >( baseIndex[idim] ) / factor) );
        }

      outputRegion.SetIndex(outputIndex);
      outputRegion.SetSize(outputSize);

      // make sure the region is within the largest possible region
      outputRegion.Crop( this->GetOutput(ilevel)->
                         GetLargestPossibleRegion() );
      // set the requested region
      this->GetOutput(ilevel)->SetRequestedRegion(outputRegion);
      }
    }
}

/**
 * GenerateInputRequestedRegion
 */
template< typename TInputImage, typename TOutputImage, typename TWaveletFunction >
void
InverseWaveletFTFilterBank< TInputImage, TOutputImage, TWaveletFunction >
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  InputImagePointer inputPtr =
    const_cast< InputImageType * >( this->GetInput() );
  if ( !inputPtr )
    {
    itkExceptionMacro(<< "Input has not been set.");
    }

  // compute baseIndex and baseSize
  typedef typename OutputImageType::SizeType   SizeType;
  typedef typename OutputImageType::IndexType  IndexType;
  typedef typename OutputImageType::RegionType RegionType;

  unsigned int refLevel = m_NumberOfLevels - 1;
  SizeType     baseSize = this->GetOutput(refLevel)->GetRequestedRegion().GetSize();
  IndexType    baseIndex = this->GetOutput(refLevel)->GetRequestedRegion().GetIndex();
  RegionType   baseRegion;

  unsigned int idim;
  for ( idim = 0; idim < ImageDimension; idim++ )
    {
    unsigned int factor = m_Schedule[refLevel][idim];
    baseIndex[idim] *= static_cast< IndexValueType >( factor );
    baseSize[idim] *= static_cast< SizeValueType >( factor );
    }
  baseRegion.SetIndex(baseIndex);
  baseRegion.SetSize(baseSize);

  // compute requirements for the smoothing part
  typedef typename TOutputImage::PixelType                    OutputPixelType;
  typedef GaussianOperator< OutputPixelType, ImageDimension > OperatorType;

  OperatorType *oper = new OperatorType;

  typename TInputImage::SizeType radius;

  RegionType inputRequestedRegion = baseRegion;
  refLevel = 0;

  for ( idim = 0; idim < TInputImage::ImageDimension; idim++ )
    {
    oper->SetDirection(idim);
    oper->SetVariance( itk::Math::sqr( 0.5 * static_cast< float >(
                                       m_Schedule[refLevel][idim] ) ) );
    oper->SetMaximumError(m_MaximumError);
    oper->CreateDirectional();
    radius[idim] = oper->GetRadius()[idim];
    }
  delete oper;

  inputRequestedRegion.PadByRadius(radius);

  // make sure the requested region is within the largest possible
  inputRequestedRegion.Crop( inputPtr->GetLargestPossibleRegion() );

  // set the input requested region
  inputPtr->SetRequestedRegion(inputRequestedRegion);
}

} // end namespace itk
#endif
