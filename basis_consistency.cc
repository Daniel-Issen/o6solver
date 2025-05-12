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

#include "basis_consistency.h"
#include "pairing.h"
#include <tuple>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>

// lookup tables to eliminate conditional logic.

// masks to update a basis (i,j,k) based on the state of the pairs
// (i,j), (i,k), (j,k).

static constexpr uint8_t ij_basis_clear_masks[16] =
  { 0,
    SET_NEG_NEG_ANY,
    SET_NEG_POS_ANY,
    SET_NEG_ANY_ANY,
    SET_POS_NEG_ANY,
    (SET_POS_NEG_ANY | SET_NEG_NEG_ANY),
    (SET_POS_NEG_ANY | SET_NEG_POS_ANY),
    (SET_POS_NEG_ANY | SET_NEG_ANY_ANY),
    SET_POS_POS_ANY,
    (SET_POS_POS_ANY | SET_NEG_NEG_ANY),
    (SET_POS_POS_ANY | SET_NEG_POS_ANY),
    (SET_POS_POS_ANY | SET_NEG_ANY_ANY),
    (SET_POS_POS_ANY | SET_POS_NEG_ANY),
    (SET_POS_POS_ANY | SET_POS_NEG_ANY | SET_NEG_NEG_ANY),
    (SET_POS_POS_ANY | SET_POS_NEG_ANY | SET_NEG_POS_ANY),
    (SET_POS_POS_ANY | SET_POS_NEG_ANY | SET_NEG_POS_ANY | SET_NEG_NEG_ANY)
  };
static constexpr uint8_t ik_basis_clear_masks[16] =
  { 0,
    SET_NEG_ANY_NEG,
    SET_NEG_ANY_POS,
    SET_NEG_ANY_ANY,
    SET_POS_ANY_NEG,
    (SET_POS_ANY_NEG | SET_NEG_ANY_NEG),
    (SET_POS_ANY_NEG | SET_NEG_ANY_POS),
    (SET_POS_ANY_NEG | SET_NEG_ANY_ANY),
    SET_POS_ANY_POS,
    (SET_POS_ANY_POS | SET_NEG_ANY_NEG),
    (SET_POS_ANY_POS | SET_NEG_ANY_POS),
    (SET_POS_ANY_POS | SET_NEG_ANY_ANY),
    (SET_POS_ANY_POS | SET_POS_ANY_NEG),
    (SET_POS_ANY_POS | SET_POS_ANY_NEG | SET_NEG_ANY_NEG),
    (SET_POS_ANY_POS | SET_POS_ANY_NEG | SET_NEG_ANY_POS),
    (SET_POS_ANY_POS | SET_POS_ANY_NEG | SET_NEG_ANY_POS | SET_NEG_ANY_NEG)
  };

static constexpr uint8_t jk_basis_clear_masks[16] =
  { 0,
    SET_ANY_NEG_NEG,
    SET_ANY_NEG_POS,
    SET_ANY_NEG_ANY,
    SET_ANY_POS_NEG,
    (SET_ANY_POS_NEG | SET_ANY_NEG_NEG),
    (SET_ANY_POS_NEG | SET_ANY_NEG_POS),
    (SET_ANY_POS_NEG | SET_ANY_NEG_ANY),
    SET_ANY_POS_POS,
    (SET_ANY_POS_POS | SET_ANY_NEG_NEG),
    (SET_ANY_POS_POS | SET_ANY_NEG_POS),
    (SET_ANY_POS_POS | SET_ANY_NEG_ANY),
    (SET_ANY_POS_POS | SET_ANY_POS_NEG),
    (SET_ANY_POS_POS | SET_ANY_POS_NEG | SET_ANY_NEG_NEG),
    (SET_ANY_POS_POS | SET_ANY_POS_NEG | SET_ANY_NEG_POS),
    (SET_ANY_POS_POS | SET_ANY_POS_NEG | SET_ANY_NEG_POS | SET_ANY_NEG_NEG)
  };

// given a basis state, bit masks to update the pairs and terms
static constexpr uint8_t basis_to_pair_and_term_clear_masks[256][6] =
  { {0, 0, 0, 0, 0, 0},
    { (SET_NEG_NEG ), 
      (SET_NEG_NEG ), 
      (SET_NEG_NEG ), 
      (SET_NEG ), 
      (SET_NEG ), 
      (SET_NEG ) },
    { (SET_NEG_NEG ), 
      (SET_NEG_POS ), 
      (SET_NEG_POS ), 
      (SET_NEG ), 
      (SET_NEG ), 
      (SET_POS ) },
    { (SET_NEG_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS ), 
      (SET_NEG_NEG ), 
      (SET_POS_NEG ), 
      (SET_NEG ), 
      (SET_POS ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_POS | SET_NEG_NEG ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS ), 
      (SET_NEG_POS ), 
      (SET_POS_POS ), 
      (SET_NEG ), 
      (SET_POS ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_POS ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_POS | SET_NEG_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_NEG ), 
      (SET_POS_NEG ), 
      (SET_NEG_NEG ), 
      (SET_POS ), 
      (SET_NEG ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_POS_NEG | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_POS_NEG | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_NEG ), 
      (SET_POS_POS ), 
      (SET_NEG_POS ), 
      (SET_POS ), 
      (SET_NEG ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_NEG ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_POS ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_POS | SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS | SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_POS ), 
      (SET_POS_NEG ), 
      (SET_POS_NEG ), 
      (SET_POS ), 
      (SET_POS ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_POS_NEG | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_POS_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_POS_NEG | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_POS_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_POS_POS | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_POS | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_POS_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_POS | SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS | SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_POS ), 
      (SET_POS_POS ), 
      (SET_POS_POS ), 
      (SET_POS ), 
      (SET_POS ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_POS | SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS | SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_POS ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_POS_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_POS_POS | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_POS | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_POS_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_NEG_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_NEG_POS | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_NEG_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_POS | SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_POS | SET_NEG_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_POS | SET_NEG_NEG | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_POS | SET_POS_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_POS_NEG | SET_POS_POS | SET_NEG_NEG | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_POS_NEG | SET_POS_POS | SET_NEG_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_NEG_NEG | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_POS | SET_POS_NEG | SET_POS_POS | SET_NEG_NEG ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_POS | SET_NEG ) },
    { (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG_NEG | SET_NEG_POS | SET_POS_NEG | SET_POS_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ), 
      (SET_NEG | SET_POS ) } };

static constexpr uint8_t threed_intermediary_set_masks[3][3][3] =
  { { { SET_NONE, SET_NONE, SET_NONE },
      { SET_NONE, SET_NONE, SET_NONE },
      { SET_NONE, SET_NONE, SET_NONE } },
    { { SET_NONE, SET_NONE, SET_NONE },
      { SET_NONE, SET_NEG_NEG_NEG, SET_NEG_NEG_POS },
      { SET_NONE, SET_NEG_POS_NEG, SET_NEG_POS_POS },},
    { { SET_NONE, SET_NONE, SET_NONE },
      { SET_NONE, SET_POS_NEG_NEG, SET_POS_NEG_POS },
      { SET_NONE, SET_POS_POS_NEG, SET_POS_POS_POS },},};

// Helper function to print a term state
std::string term_state_str(uint8_t state) {
  switch (state) {
  case 0: return "CONTRADICTION";
  case SET_NEG: return "NEG (0x1)";
  case SET_POS: return "POS (0x2)";
  case SET_ANY: return "ANY (0x3)";
  default: return "UNKNOWN (" + std::to_string(state) + ")";
  }
}

// Helper function to print a pair state
std::string pair_state_str(uint8_t state) {
  std::string result = "0x" + std::to_string(state) + " [";
  if (state & SET_NEG_NEG) result += "NEG-NEG ";
  if (state & SET_NEG_POS) result += "NEG-POS ";
  if (state & SET_POS_NEG) result += "POS-NEG ";
  if (state & SET_POS_POS) result += "POS-POS";
  result += "]";
  return result;
}

// Helper function to print a basis state
std::string basis_state_str(uint8_t state) {
  if (!state) return "CONTRADICTION";
    
  std::string result = "0x" + std::to_string(state) + " [";
  if (state & SET_NEG_NEG_NEG) result += "NEG-NEG-NEG ";
  if (state & SET_NEG_NEG_POS) result += "NEG-NEG-POS ";
  if (state & SET_NEG_POS_NEG) result += "NEG-POS-NEG ";
  if (state & SET_NEG_POS_POS) result += "NEG-POS-POS ";
  if (state & SET_POS_NEG_NEG) result += "POS-NEG-NEG ";
  if (state & SET_POS_NEG_POS) result += "POS-NEG-POS ";
  if (state & SET_POS_POS_NEG) result += "POS-POS-NEG ";
  if (state & SET_POS_POS_POS) result += "POS-POS-POS";
  result += "]";
  return result;
}

// propagate term states to pairs then propage the pair
// state back down to the terms.
UpdateResult update_pair_states(uint64_t i, uint64_t j,
				std::vector<uint8_t>& term_states,
				std::vector<uint8_t>& pair_states) {
  // Save original states to detect changes
  uint8_t &term_i = term_states[i];
  uint8_t term_i_orig = term_i;

  uint8_t &term_j = term_states[j];
  uint8_t term_j_orig = term_j;

  // Get pair indices
  uint64_t ij_idx = pair2d(i, j);
    
  // Save original pair states
  uint8_t &pair_ij = pair_states[ij_idx];
  uint8_t pair_ij_orig = pair_ij;

  // update the pairs and basis from the terms
  switch (term_i) {
  case 0:
    // UNSAT
    return UpdateResult(true,true);
    break;
  case CLEAR_POS:
    pair_ij &= CLEAR_POS_ANY;
    break;
  case CLEAR_NEG:
    pair_ij &= CLEAR_NEG_ANY;
    break;
  }

  switch(term_j) {
  case 0:
    // UNSAT
    return UpdateResult(true,true);
    break;
  case CLEAR_POS:
    pair_ij &= CLEAR_ANY_POS;
    break;
  case CLEAR_NEG:
    pair_ij &= CLEAR_ANY_NEG;
    break;
  }
  if(!pair_ij) {
    // UNSAT
    return UpdateResult(true,true);
  }
  if(!(pair_ij & SET_NEG_ANY)) {
    term_i &= CLEAR_NEG;
  }
  if(!(pair_ij & SET_POS_ANY)) {
    term_i &= CLEAR_POS;
  }
  if(!term_i) {
    return UpdateResult(true,true);
  }
  if(!(pair_ij & SET_ANY_NEG)) {
    term_j &= CLEAR_NEG;
  }
  if(!(pair_ij & SET_ANY_POS)) {
    term_j &= CLEAR_POS;
  }
  if(!term_j) {
    return UpdateResult(true,true);
  }

  if((pair_ij != pair_ij_orig)    ||
     (term_i  != term_i_orig)     ||
     (term_j  != term_j_orig)) {
    return UpdateResult(true, false);
  }
  return UpdateResult(false, false);
}

// propagate term states to pairs to a basis then propage the basis
// state back down to the pairs and terms.
UpdateResult update_basis_states(uint64_t i, uint64_t j, uint64_t k,
				 uint64_t basis_idx,
				 std::vector<uint8_t>& term_states,
				 std::vector<uint8_t>& pair_states,
				 std::vector<uint8_t>& basis_states) {
  // Save original states to detect changes
  uint8_t &term_i = term_states[i];
  uint8_t term_i_orig = term_i;

  uint8_t &term_j = term_states[j];
  uint8_t term_j_orig = term_j;

  uint8_t &term_k = term_states[k];
  uint8_t term_k_orig = term_k;
    
  // Get pair indices
  uint64_t ij_idx = pair2d(i, j);
  uint64_t ik_idx = pair2d(i, k);
  uint64_t jk_idx = pair2d(j, k);
    
  // Save original pair states
  uint8_t &pair_ij = pair_states[ij_idx];
  uint8_t pair_ij_orig = pair_ij;

  uint8_t &pair_ik = pair_states[ik_idx];
  uint8_t pair_ik_orig = pair_ik;
  
  uint8_t &pair_jk = pair_states[jk_idx];
  uint8_t pair_jk_orig = pair_jk;
    
  // Save original basis state
  uint8_t &basis_ijk = basis_states[basis_idx];
  uint8_t basis_ijk_orig = basis_ijk;

  // update the pairs and basis from the terms
  switch (term_i) {
  case 0:
    // UNSAT
    return UpdateResult(true,true);
    break;
  case CLEAR_POS:
    pair_ij &= CLEAR_POS_ANY;
    pair_ik &= CLEAR_POS_ANY;
    basis_ijk &= CLEAR_POS_ANY_ANY;
    break;
  case CLEAR_NEG:
    pair_ij &= CLEAR_NEG_ANY;
    pair_ik &= CLEAR_NEG_ANY;
    basis_ijk &= CLEAR_NEG_ANY_ANY;
    break;
  }

  switch(term_j) {
  case 0:
    // UNSAT
    return UpdateResult(true,true);
    break;
  case CLEAR_POS:
    pair_ij &= CLEAR_ANY_POS;
    pair_jk &= CLEAR_POS_ANY;
    basis_ijk &= CLEAR_ANY_POS_ANY;
    break;
  case CLEAR_NEG:
    pair_ij &= CLEAR_ANY_NEG;
    pair_jk &= CLEAR_NEG_ANY;
    basis_ijk &= CLEAR_ANY_NEG_ANY;
    break;
  }
  
  switch(term_k) {
  case 0:
    // UNSAT
    return UpdateResult(true,true);
    break;
  case CLEAR_POS:
    pair_ik &= CLEAR_ANY_POS;
    pair_jk &= CLEAR_ANY_POS;
    basis_ijk &= CLEAR_ANY_ANY_POS;
    break;
  case CLEAR_NEG:
    pair_ik &= CLEAR_ANY_NEG;
    pair_jk &= CLEAR_ANY_NEG;
    basis_ijk &= CLEAR_ANY_ANY_NEG;
    break;
  }
  
  // update basis from the pairs
  basis_ijk &= ij_basis_clear_masks[pair_ij];
  basis_ijk &= ik_basis_clear_masks[pair_ik];
  basis_ijk &= jk_basis_clear_masks[pair_jk];

  // update pairs and terms from the basis
  const uint8_t *bpt_clear_masks =
    basis_to_pair_and_term_clear_masks[basis_ijk];
  pair_ij &= bpt_clear_masks[0];
  pair_ik &= bpt_clear_masks[1];
  pair_jk &= bpt_clear_masks[2];
  term_i  &= bpt_clear_masks[3];
  term_j  &= bpt_clear_masks[4];
  term_k  &= bpt_clear_masks[5];

  if((!basis_ijk) ||
     (!pair_ij)   ||
     (!pair_ik)   ||
     (!pair_jk)   ||
     (!term_i)    ||
     (!term_j)    ||
     (!term_k)) {
    // UNSAT
    return UpdateResult(true,true);
  }
  if((basis_ijk != basis_ijk_orig)||
     (pair_ij != pair_ij_orig)    ||
     (pair_ik != pair_ik_orig)    ||
     (pair_jk != pair_jk_orig)    ||
     (term_i  != term_i_orig)     ||
     (term_j  != term_j_orig)     ||
     (term_k  != term_k_orig)) {
    return UpdateResult(true, false);
  }
  return UpdateResult(false, false);
}

// First, define a fixed-size container for intermediary bases
struct IntermediaryBasis {
  uint64_t basis_idx;
  uint64_t i;
  uint64_t j;
  uint64_t k;
  unsigned int offset1, offset2, offset3;
  uint8_t state;
};

// Maximum number of intermediary bases: (6 choose 3) - 2 original
// bases = 18
constexpr size_t MAX_INTERMEDIARY_BASES = 18;

// Intermediary is a tuple of the basis index, the three term
// offsets, and a new basis state.  intermediary bases are composed
// of one term from one basis and two terms from the other.  the key
// insight here is that these intermediaries are first made
// conconsistent and representative (wrt consistency with basis1 and
// basis2) of all other bases that share terms by calling
// update_basis_states.  A conflict with a representative
// intermediary basis (a,b,c) is representative of all possible
// other bases (a,b,X), (a,X,c), (X,b,c),(a,X,Y),(X,b,Y),(X,Y,c) and
// once basis(a,b,c) is trimmed of states that are inconsistent with
// basis1, basis2, and all the other representative intermediary
// bases, the global state is made consistent by calling
// update_basis_states on (a,b,c).
void generate_intermediaries(const uint64_t* b1_array,
			     const uint64_t* b2_array,
			     IntermediaryBasis* intermediaries,
			     size_t* num_intermediaries) {
    
  // Merged array of terms and their offsets
  uint64_t all_terms[6];
  unsigned int all_offsets[6];
  size_t term_count = 0;
    
  // Parallel traversal of both sorted arrays (same as before)
  int i = 0, j = 0;
  while (i < 3 && j < 3) {
    if (b1_array[i] < b2_array[j]) {
      all_terms[term_count] = b1_array[i];
      all_offsets[term_count] = i;
      term_count++;
      i++;
    } else if (b1_array[i] > b2_array[j]) {
      all_terms[term_count] = b2_array[j];
      all_offsets[term_count] = j + 3;
      term_count++;
      j++;
    } else {
      all_terms[term_count] = b1_array[i];
      all_offsets[term_count] = i;
      term_count++;
      i++;
      j++;
    }
  }
    
  // Add any remaining terms
  while (i < 3) {
    all_terms[term_count] = b1_array[i];
    all_offsets[term_count] = i;
    term_count++;
    i++;
  }
    
  while (j < 3) {
    all_terms[term_count] = b2_array[j];
    all_offsets[term_count] = j + 3;
    term_count++;
    j++;
  }
    
  // Create indices for basis1 and basis2 to filter them out later
  uint64_t basis1_idx = pair3d(b1_array[0], b1_array[1], b1_array[2]);
  uint64_t basis2_idx = pair3d(b2_array[0], b2_array[1], b2_array[2]);
    
  // Reset counter
  *num_intermediaries = 0;
    
  // Generate all combinations of 3 terms from the merged list
  for (size_t i = 0;
       ((i < term_count) &&
	(*num_intermediaries < MAX_INTERMEDIARY_BASES));
       i++) {
    for (size_t j = i + 1;
	 ((j < term_count) &&
	  (*num_intermediaries < MAX_INTERMEDIARY_BASES));
	 j++) {
      for (size_t k = j + 1;
	   ((k < term_count) &&
	    (*num_intermediaries < MAX_INTERMEDIARY_BASES));
	   k++) {
	// Create basis index
	uint64_t basis_idx =
	  pair3d(all_terms[i], all_terms[j], all_terms[k]);
                
	// Skip if this is one of the original bases
	if (basis_idx == basis1_idx || basis_idx == basis2_idx) {
	  continue;
	}
                
	// Add intermediary with pre-calculated offsets (using fixed
	// array)
	intermediaries[*num_intermediaries].basis_idx = basis_idx;
	intermediaries[*num_intermediaries].i = all_terms[i];
	intermediaries[*num_intermediaries].j = all_terms[j];
	intermediaries[*num_intermediaries].k = all_terms[k];
	intermediaries[*num_intermediaries].offset1 = all_offsets[i];
	intermediaries[*num_intermediaries].offset2 = all_offsets[j];
	intermediaries[*num_intermediaries].offset3 = all_offsets[k];
	intermediaries[*num_intermediaries].state = 0;
                
	(*num_intermediaries)++;
      }
    }
  }
}

UpdateResult ensure_basis_consistency
(uint64_t i1,
 uint64_t j1,
 uint64_t k1,
 uint64_t i2,
 uint64_t j2,
 uint64_t k2,
 uint64_t basis1_idx,
 uint64_t basis2_idx,
 std::vector<uint8_t>& term_states,
 std::vector<uint8_t>& pair_states,
 std::vector<uint8_t>& basis_states) {

  // First update each basis individually
  UpdateResult result = update_basis_states(i1, j1, k1,
					    basis1_idx,
					    term_states,
					    pair_states,
					    basis_states);
  
  if(result.has_zero) {
    return result;
  }
  // make basis2 consistent with basis1
  UpdateResult basis2_result = update_basis_states(i2, j2, k2,
						   basis2_idx,
						   term_states,
						   pair_states,
						   basis_states);
  if(basis2_result.has_zero) {
    return basis2_result;
  }
  result.changed = result.changed || basis2_result.changed;

  // make basis1 consistent with basis2
  result = update_basis_states(i1, j1, k1,
			       basis1_idx,
			       term_states,
			       pair_states,
			       basis_states);
  if(result.has_zero) {
    return result;
  }
    
  // Setup for intermediary generation
  uint64_t b1_array[] = {i1, j1, k1};
  uint64_t b2_array[] = {i2, j2, k2};
    
  // Use stack-allocated array for intermediaries
  IntermediaryBasis intermediaries[MAX_INTERMEDIARY_BASES];
  size_t num_intermediaries = 0;
    
  // Generate intermediaries directly from the variable triplets
  generate_intermediaries(b1_array,
			  b2_array,
			  intermediaries,
			  &num_intermediaries);
    
  // Update all intermediary bases
  bool any_changed;
  do {
    any_changed = false;
    for (size_t idx = 0; idx < num_intermediaries; idx++) {
      // We still need to unpack this one intermediary basis
      UpdateResult inter_result =
	update_basis_states(intermediaries[idx].i,
			    intermediaries[idx].j,
			    intermediaries[idx].k,
			    intermediaries[idx].basis_idx,
			    term_states,
			    pair_states,
			    basis_states);
      if (inter_result.has_zero) {
	return inter_result;
      }
      any_changed |= inter_result.changed;
    }
  } while (any_changed);
  // Calculate consistent states
  uint8_t &basis1_state = basis_states[basis1_idx];
  uint8_t &basis2_state = basis_states[basis2_idx];
    
  uint8_t new_basis1_state = 0;
  uint8_t new_basis2_state = 0;
    
  // Pre-compute pair indices to avoid recalculation
  uint64_t ij1_idx = pair2d(i1, j1);
  uint64_t ik1_idx = pair2d(i1, k1);
  uint64_t jk1_idx = pair2d(j1, k1);
  uint64_t ij2_idx = pair2d(i2, j2);
  uint64_t ik2_idx = pair2d(i2, k2);
  uint64_t jk2_idx = pair2d(j2, k2);

  // Fixed-size array for intermediary proposals
  uint8_t intermediary_proposals[MAX_INTERMEDIARY_BASES] = {0};
    
  // For each set bit in basis1_state
  uint8_t basis1_bits = basis1_state;
  while (basis1_bits) {
    uint8_t basis1_bit =
      basis1_bits & -basis1_bits;  // Extract lowest set bit
    basis1_bits &= ~basis1_bit;	   // Clear that bit
        
    const uint8_t *first_bpt_clear_masks =
      basis_to_pair_and_term_clear_masks[basis1_bit];
        
    // For each set bit in basis2_state
    uint8_t basis2_bits = basis2_state;
    while (basis2_bits) {
      uint8_t basis2_bit = basis2_bits & -basis2_bits;
      basis2_bits &= ~basis2_bit;
                    
      const uint8_t *second_bpt_clear_masks =
	basis_to_pair_and_term_clear_masks[basis2_bit];
                    
      // Get joint term states for this pair of basis states
      uint8_t joint_states[6] = {
	first_bpt_clear_masks[3],
	first_bpt_clear_masks[4],
	first_bpt_clear_masks[5],
	second_bpt_clear_masks[3],
	second_bpt_clear_masks[4],
	second_bpt_clear_masks[5],
      };
                    
      // Calculate all required states first for better memory
      // locality
      for (size_t i = 0; i < num_intermediaries; i++) {
	uint8_t i_state = joint_states[intermediaries[i].offset1];
	uint8_t j_state = joint_states[intermediaries[i].offset2];
	uint8_t k_state = joint_states[intermediaries[i].offset3];
  
	intermediary_proposals[i] =
	  threed_intermediary_set_masks[i_state][j_state][k_state];
      }
      
      // Check if this combination is consistent with all
      // intermediaries
      bool consistent = true;
      for (size_t i = 0; i < num_intermediaries; i++) {
	uint64_t inter_idx = intermediaries[i].basis_idx;
	if (!(basis_states[inter_idx] & intermediary_proposals[i])) {
	  consistent = false;
	  break;
	}
      }
                    
      if (consistent) {
	// This pair of basis states is consistent
	new_basis1_state |= basis1_bit;
	new_basis2_state |= basis2_bit;
                        
	// Update all intermediary state values
	for (size_t i = 0; i < num_intermediaries; i++) {
	  intermediaries[i].state |= intermediary_proposals[i];
	}
      }
    }
  }

  // Update basis1 if changed
  if(basis1_state != new_basis1_state) {
    basis1_state = new_basis1_state;
    result.changed = true;
        
    // Update pairs and terms using pre-computed indices
    const uint8_t *bpt_clear_masks =
      basis_to_pair_and_term_clear_masks[new_basis1_state];

    pair_states[ij1_idx] &= bpt_clear_masks[0];
    pair_states[ik1_idx] &= bpt_clear_masks[1];
    pair_states[jk1_idx] &= bpt_clear_masks[2];
    term_states[i1] &= bpt_clear_masks[3];
    term_states[j1] &= bpt_clear_masks[4];
    term_states[k1] &= bpt_clear_masks[5];

    // Check for contradictions
    if((!new_basis1_state) ||
       (!pair_states[ij1_idx]) ||
       (!pair_states[ik1_idx]) ||
       (!pair_states[jk1_idx]) ||
       (!term_states[i1]) ||
       (!term_states[j1]) ||
       (!term_states[k1])) {
      result.has_zero = true;
      return result;
    }
  }

  // Update basis2 if changed
  if(basis2_state != new_basis2_state) {
    basis2_state = new_basis2_state;
    result.changed = true;
        
    const uint8_t *bpt_clear_masks =
      basis_to_pair_and_term_clear_masks[new_basis2_state];

    pair_states[ij2_idx] &= bpt_clear_masks[0];
    pair_states[ik2_idx] &= bpt_clear_masks[1];
    pair_states[jk2_idx] &= bpt_clear_masks[2];
    term_states[i2] &= bpt_clear_masks[3];
    term_states[j2] &= bpt_clear_masks[4];
    term_states[k2] &= bpt_clear_masks[5];

    // Check for contradictions
    if((!new_basis2_state) ||
       (!pair_states[ij2_idx]) ||
       (!pair_states[ik2_idx]) ||
       (!pair_states[jk2_idx]) ||
       (!term_states[i2]) ||
       (!term_states[j2]) ||
       (!term_states[k2])) {
      result.has_zero = true;
      return result;
    }
  }

  // Update intermediary basis states
  for(size_t i = 0; i < num_intermediaries; i++) {
    basis_states[intermediaries[i].basis_idx] = intermediaries[i].state;
  }
    
  return result;
}

void ensure_global_consistency(std::vector<uint8_t>& term_states,
			       std::vector<uint8_t>& pair_states,
			       std::vector<uint8_t>& basis_states,
			       bool& has_contradiction,
			       uint64_t starting_basis_pair,
			       uint64_t ending_basis_pair) {
  has_contradiction = false;
  bool changed = true;

  while (changed) {
    changed = false;
    for(uint64_t basis_pair = starting_basis_pair;
	basis_pair < ending_basis_pair;
	++basis_pair) {
      uint64_t basis1_idx,basis2_idx;
      std::tie(basis1_idx,basis2_idx) = unpair2d(basis_pair);
      uint64_t i1,j1,k1,i2,j2,k2;
      std::tie(i1,j1,k1) = unpair3d(basis1_idx);
      std::tie(i2,j2,k2) = unpair3d(basis2_idx);
      auto result =
	ensure_basis_consistency(i1, j1, k1, i2, j2, k2,
				 basis1_idx,
				 basis2_idx,
				 term_states,
				 pair_states,
				 basis_states);
      if (result.has_zero) {
	has_contradiction = true;
	return;
      }
                        
      if (result.changed) {
	changed = true;
      }
    }
  }
  return;
}  
