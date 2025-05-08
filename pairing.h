// MIT License

// Copyright (c) 2025 Daniel Issen

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// pairing 2d and 3d functions.  see
// https://en.wikipedia.org/wiki/Pairing_function.  for our problem we
// are only interested in the unique combinations of terms where no
// terms are repeated (i < j < k) so we fill arrays in upper triagular
// form and avoid wasting space for the redundant combinations.

#pragma once
#include <tuple>
#include <cstdint>

// Maps (i,j) where i < j to a unique index
uint64_t pair2d(uint64_t i, uint64_t j);

// Maps index back to (i,j) where i < j
std::tuple<uint64_t, uint64_t> unpair2d(uint64_t index);

// Maps (i,j,k) where i < j < k to a unique index
uint64_t pair3d(uint64_t i, uint64_t j, uint64_t k);

// Maps index back to (i,j,k) where i < j < k
std::tuple<uint64_t, uint64_t, uint64_t> unpair3d(uint64_t index);

// Helper functions to calculate total array sizes
uint64_t calculate_array_size_2d(uint64_t n); // N Choose 2
uint64_t calculate_array_size_3d(uint64_t n); // N Choose 3
