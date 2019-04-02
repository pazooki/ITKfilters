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

#include <itkImage.h>
#include <itkChangeInformationImageFilter.h>

// do not forget to call Update() in the filter generating imageToRestore
template<typename TImage, typename TReferenceImage>
void restoreMetadata(TImage * imageToRestore, const TReferenceImage * referenceImage, size_t level)
{
  auto spacing = referenceImage->GetSpacing();
  auto origin = referenceImage->GetOrigin();
  auto direction = referenceImage->GetDirection();
  imageToRestore->SetSpacing(spacing * std::pow(2, level) );
  imageToRestore->SetOrigin(origin);
  imageToRestore->SetDirection(direction);
}

template<typename TImage, typename TReferenceImage>
typename TImage::Pointer changeMetadata(TImage* imageToRestore, const TReferenceImage * referenceImage, size_t level)
{
  using ChangeInfoFilter = itk::ChangeInformationImageFilter<TImage>;
  auto spacing = referenceImage->GetSpacing();
  auto origin = referenceImage->GetOrigin();
  auto direction = referenceImage->GetDirection();
  std::cout << "Reference Image: " << std::endl;
  std::cout << origin << std::endl;
  std::cout << spacing << std::endl;
  std::cout << direction << std::endl;
  auto changeFilter = ChangeInfoFilter::New();
  changeFilter->SetInput(imageToRestore);
  changeFilter->ChangeDirectionOn();
  changeFilter->ChangeOriginOn();
  changeFilter->ChangeSpacingOn();
  changeFilter->SetOutputOrigin(origin);
  changeFilter->SetOutputSpacing(spacing * std::pow(2, level) );
  changeFilter->SetOutputDirection(direction);
  changeFilter->Update();
  return changeFilter->GetOutput();
}

