/*=========================================================================
 *  Linear index to subindices.
 *  Copyright 2016 Pablo Hernandez-Cerdan
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

#ifndef Ind2Sub_h
#define Ind2Sub_h
#include <array>
#include <iostream>

/**
 * @brief Return SubIndex [i,j,k,...] from linear index and matrix size.
 * Based on Ind2Sub from Matlab.
 * @note Indices range from 0 to ns - 1
 * Linear index cannot exceed the accumulated product of sizes, throwing runtime errors.
 *
 * @tparam D Dimension of the matrix
 * @param linear_index linear index, for example in 3D linear_index = i + nx(j + ny*k)
 * @param ns array containing sizes per dimenstion (nx,ny,...)
 *
 * @return array with subindexes: [i,j,k, ...
 */
template <unsigned int D>
std::array<unsigned int, D>
Ind2Sub(const unsigned int &linear_index, const std::array<unsigned int, D> & ns)
{
  // accumulative product.
  std::array<unsigned int, D> cumprod;
  unsigned int accum = 1;
  cumprod[0] = accum;
  for (unsigned int i = 1 ; i < D  ; ++i)
    {
    accum *= ns[i -1] ;
    cumprod[i] = accum;
    }
  unsigned int max_index = accum * ns[D -1] - 1;
  if (linear_index > max_index)
    throw std::runtime_error("Ind2Sub: input index is incompatible with the given size");

  std::array<unsigned int, D> out ;
  unsigned int i(0), rem(0), temp_index(linear_index);
  for (int si = D - 1 ; si > -1 ; --si)
    {
    i = si;
    rem = (temp_index ) % cumprod[i] ;
    out[i] = (temp_index - rem) / cumprod[i] ;
    temp_index = rem;
    }

  return out;

};
#endif

// //TEST
// int main(int, char**)
// {
//   const unsigned int D = 3;
//   std::array<unsigned int, D> sizes{ {2,2,2} };
//   std::cout << "sizes: (2,2,2) "<< std::endl;
//   for (unsigned int n = 0 ; n < 2*2*2 ; ++n)
//     {
//     std::cout << "index " << n  << " -> ";
//     auto out = Ind2Sub<D>(n, sizes);
//     for (unsigned int i = 0 ; i < D ; ++i)
//       std::cout <<out[i] << " , ";
//     std::cout << '\n';
//
//     }
//
//   const unsigned int D2 = 2;
//   std::array<unsigned int, D2> sizes2{ {3,4} };
//   std::cout << "sizes: (3,4) "<< std::endl;
//   for (unsigned int n = 0 ; n < 3*4 ; ++n)
//     {
//     std::cout << "index " << n  << " -> ";
//     auto out = Ind2Sub<D2>(n, sizes2);
//     for (unsigned int i = 0 ; i < D2 ; ++i)
//       std::cout <<out[i] << " , ";
//     std::cout << '\n';
//
//     }
// }
