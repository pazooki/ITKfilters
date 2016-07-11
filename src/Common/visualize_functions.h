#pragma once
#include <cstddef>
namespace visualize
{
    template<typename T >
    void VisualizeITKImage(const T* img, size_t win_x = 600, size_t win_y = 600);

    template<typename TLeft, typename TRight >
    void VisualizeITKImages(const TLeft* leftImg, const TRight* rightImg,
                              size_t win_x = 800, size_t win_y = 800);
}// visualize
#include "visualize_functions.hxx"
