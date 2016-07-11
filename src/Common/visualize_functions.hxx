#pragma once
#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleRubberBand3D.h>
#include <vtkRenderer.h>
#include <vtkImageMapper.h>
#include <vtkImagePlaneWidget.h>
#include "itkImage.h"
#include "itkImageToVTKImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include <array>
namespace visualize
{
template<typename T >
void VisualizeITKImage(const T* img, size_t win_x , size_t win_y )
{
      const auto kDimension = T::ImageDimension;
      using ConnectorType = itk::ImageToVTKImageFilter<T> ;
      auto connector = ConnectorType::New();
      connector->SetInput(img);
      connector->Update();
      connector->UpdateLargestPossibleRegion();

      // Setup renderers
      auto renderer = vtkSmartPointer<vtkRenderer>::New();

      // Setup render window
      auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
      renderWindow->SetSize(win_x, win_y);
      renderWindow->AddRenderer(renderer);

      // Setup render window interactor
      auto renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
      auto style = vtkSmartPointer<vtkInteractorStyleRubberBand3D>::New();
      renderWindowInteractor->SetInteractorStyle(style);

      // Render and start interaction
      renderWindowInteractor->SetRenderWindow(renderWindow);

      //Prepare for slices.
      using FilterType = itk::StatisticsImageFilter<T> ;
      auto filter = FilterType::New();
      filter->SetInput(img);
      filter->Update();
      filter->UpdateLargestPossibleRegion();
      double min_intensity = filter->GetMinimum();
      double max_intensity = filter->GetMaximum();
      auto window = max_intensity - min_intensity;
      auto level = min_intensity + window / 2;
    /** SLICES */
      std::array<vtkSmartPointer<vtkImagePlaneWidget>, kDimension> slice_planes;
      for (unsigned i = 0; i < kDimension; ++i) {
          slice_planes[i] = vtkSmartPointer<vtkImagePlaneWidget>::New();
          slice_planes[i]->SetResliceInterpolateToCubic();
          slice_planes[i]->DisplayTextOn();
          slice_planes[i]->SetInteractor(renderWindowInteractor);
          slice_planes[i]->PlaceWidget();
          slice_planes[i]->SetSliceIndex(0);
          slice_planes[i]->SetMarginSizeX(0);
          slice_planes[i]->SetMarginSizeY(0);
          slice_planes[i]->SetRightButtonAction(
                  vtkImagePlaneWidget::VTK_SLICE_MOTION_ACTION);
          slice_planes[i]->SetMiddleButtonAction(
                  vtkImagePlaneWidget::VTK_WINDOW_LEVEL_ACTION);
          slice_planes[i]->TextureInterpolateOff();

          slice_planes[i]->SetInputData(connector->GetOutput());
          slice_planes[i]->SetPlaneOrientation(i);
          slice_planes[i]->UpdatePlacement();
          slice_planes[i]->SetWindowLevel(window, level);
          slice_planes[i]->On();
      }
      renderer->ResetCamera();

      renderWindowInteractor->Initialize();
      renderWindowInteractor->Start();

};

template<typename TLeft, typename TRight >
void VisualizeITKImages(const TLeft* leftImg, const TRight* rightImg,
            size_t win_x , size_t win_y ){
      const auto leftDimension = TLeft::ImageDimension;
      const auto rightDimension = TRight::ImageDimension;


      using LeftConnectorType = itk::ImageToVTKImageFilter<TLeft> ;
      auto leftConnector = LeftConnectorType::New();
      leftConnector->SetInput(leftImg);
      leftConnector->Update();
      leftConnector->UpdateLargestPossibleRegion();

      using RightConnectorType = itk::ImageToVTKImageFilter<TRight> ;
      auto rightConnector = RightConnectorType::New();
      rightConnector->SetInput(rightImg);
      rightConnector->Update();
      rightConnector->UpdateLargestPossibleRegion();

      // Setup renderer (UNIQUE)
      auto renderer  = vtkSmartPointer<vtkRenderer>::New();

      // Setup render window (UNIQUE)
      auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
      renderWindow->SetSize(win_x, win_y);
      renderWindow->AddRenderer(renderer);

      // Setup render window interactor (UNIQUE)
      auto renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
      auto style = vtkSmartPointer<vtkInteractorStyleRubberBand3D>::New();
      renderWindowInteractor->SetInteractorStyle(style);

      // Render and start interaction
      renderWindowInteractor->SetRenderWindow(renderWindow);

      using LeftFilterType = itk::StatisticsImageFilter<TLeft> ;
      //Prepare for slices (BOTH)
      auto leftFilter = LeftFilterType::New();
      leftFilter->SetInput(leftImg);
      leftFilter->Update();
      leftFilter->UpdateLargestPossibleRegion();
      double leftMin_intensity = leftFilter->GetMinimum();
      double leftMax_intensity = leftFilter->GetMaximum();
      auto leftWindow = leftMax_intensity - leftMin_intensity;
      auto leftLevel = leftMin_intensity + leftWindow / 2;

      using RightFilterType = itk::StatisticsImageFilter<TRight> ;
      auto rightFilter = RightFilterType::New();
      rightFilter->SetInput(rightImg);
      rightFilter->Update();
      rightFilter->UpdateLargestPossibleRegion();
      double rightMin_intensity = rightFilter->GetMinimum();
      double rightMax_intensity = rightFilter->GetMaximum();
      auto rightWindow = rightMax_intensity - rightMin_intensity;
      auto rightLevel = rightMin_intensity + rightWindow / 2;
      /** SLICES (BOTH) */
      std::array<vtkSmartPointer<vtkImagePlaneWidget>, leftDimension> leftSlice_planes;
      for (unsigned i = 0; i < leftDimension; ++i) {
          leftSlice_planes[i] = vtkSmartPointer<vtkImagePlaneWidget>::New();
          leftSlice_planes[i]->SetResliceInterpolateToCubic();
          leftSlice_planes[i]->DisplayTextOn();
          leftSlice_planes[i]->SetInteractor(renderWindowInteractor);
          leftSlice_planes[i]->PlaceWidget();
          leftSlice_planes[i]->SetSliceIndex(0);
          leftSlice_planes[i]->SetMarginSizeX(0);
          leftSlice_planes[i]->SetMarginSizeY(0);
          leftSlice_planes[i]->SetRightButtonAction(
                  vtkImagePlaneWidget::VTK_SLICE_MOTION_ACTION);
          leftSlice_planes[i]->SetMiddleButtonAction(
                  vtkImagePlaneWidget::VTK_WINDOW_LEVEL_ACTION);
          leftSlice_planes[i]->TextureInterpolateOff();

          leftSlice_planes[i]->SetInputData(leftConnector->GetOutput());
          leftSlice_planes[i]->SetPlaneOrientation(i);
          leftSlice_planes[i]->UpdatePlacement();
          leftSlice_planes[i]->SetWindowLevel(leftWindow, leftLevel);
          leftSlice_planes[i]->On();
      }
      std::array<vtkSmartPointer<vtkImagePlaneWidget>, rightDimension> rightSlice_planes;
      for (unsigned i = 0; i < rightDimension; ++i) {
          rightSlice_planes[i] = vtkSmartPointer<vtkImagePlaneWidget>::New();
          rightSlice_planes[i]->SetResliceInterpolateToCubic();
          rightSlice_planes[i]->DisplayTextOn();
          rightSlice_planes[i]->SetInteractor(renderWindowInteractor);
          rightSlice_planes[i]->PlaceWidget();
          rightSlice_planes[i]->SetSliceIndex(0);
          rightSlice_planes[i]->SetMarginSizeX(0);
          rightSlice_planes[i]->SetMarginSizeY(0);
          rightSlice_planes[i]->SetRightButtonAction(
                  vtkImagePlaneWidget::VTK_SLICE_MOTION_ACTION);
          rightSlice_planes[i]->SetMiddleButtonAction(
                  vtkImagePlaneWidget::VTK_WINDOW_LEVEL_ACTION);
          rightSlice_planes[i]->TextureInterpolateOff();

          rightSlice_planes[i]->SetInputData(rightConnector->GetOutput());
          rightSlice_planes[i]->SetPlaneOrientation(i);
          rightSlice_planes[i]->UpdatePlacement();
          rightSlice_planes[i]->SetWindowLevel(rightWindow, rightLevel);
          rightSlice_planes[i]->On();
      }
      renderer->ResetCamera();

      renderWindowInteractor->Initialize();
      renderWindowInteractor->Start();

};
}// visualize
// #endif
