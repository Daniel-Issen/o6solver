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

// useful contants

#pragma once
#include <cstdint>

// Term state flags (2 bits - values 0-3)
static constexpr uint8_t SET_NONE = 0;
static constexpr uint8_t SET_NEG = 1;
static constexpr uint8_t SET_POS = 2;
static constexpr uint8_t SET_ANY = SET_NEG | SET_POS; // 3 - can be
						      // either NEG or
						      // POS 
static constexpr uint8_t CLEAR_NEG = 2;
static constexpr uint8_t CLEAR_POS = 1;

static constexpr uint8_t oned_clear_masks[2] = { CLEAR_NEG, CLEAR_POS };

// Pair state flags (4 bits - values 0-15)
static constexpr uint8_t PAIR_NONE = 0;
static constexpr uint8_t SET_NEG_NEG = 1;
static constexpr uint8_t SET_NEG_POS = 2;
static constexpr uint8_t SET_POS_NEG = 4;
static constexpr uint8_t SET_POS_POS = 8;
static constexpr uint8_t CLEAR_NEG_NEG = 14;
static constexpr uint8_t CLEAR_NEG_POS = 13;
static constexpr uint8_t CLEAR_POS_NEG = 11;
static constexpr uint8_t CLEAR_POS_POS = 7;

// Aggregate pair state constants 
static constexpr uint8_t SET_NEG_ANY = SET_NEG_NEG | SET_NEG_POS;     // 3 - first term NEG, second any
static constexpr uint8_t SET_POS_ANY = SET_POS_NEG | SET_POS_POS;     // 12 - first term POS, second any
static constexpr uint8_t SET_ANY_NEG = SET_NEG_NEG | SET_POS_NEG;     // 5 - first term any, second NEG
static constexpr uint8_t SET_ANY_POS = SET_NEG_POS | SET_POS_POS;     // 10 - first term any, second POS
static constexpr uint8_t SET_ANY_ANY = SET_NEG_ANY | SET_POS_ANY;     // 15 - both terms can be anything
static constexpr uint8_t CLEAR_NEG_ANY = (CLEAR_NEG_NEG & CLEAR_NEG_POS);
static constexpr uint8_t CLEAR_POS_ANY = (CLEAR_POS_NEG & CLEAR_POS_POS);
static constexpr uint8_t CLEAR_ANY_NEG = (CLEAR_NEG_NEG & CLEAR_POS_NEG);
static constexpr uint8_t CLEAR_ANY_POS = (CLEAR_NEG_POS & CLEAR_POS_POS);

static constexpr uint8_t twod_clear_masks[2][2] = {
  { CLEAR_NEG_NEG, CLEAR_NEG_POS },{ CLEAR_POS_NEG, CLEAR_POS_POS },};

// Basis state flags (8 bits - values 0-255)
static constexpr uint8_t BASIS_NONE = 0;
static constexpr uint8_t SET_NEG_NEG_NEG = 1;
static constexpr uint8_t SET_NEG_NEG_POS = 2;
static constexpr uint8_t SET_NEG_POS_NEG = 4;
static constexpr uint8_t SET_NEG_POS_POS = 8;
static constexpr uint8_t SET_POS_NEG_NEG = 16;
static constexpr uint8_t SET_POS_NEG_POS = 32;
static constexpr uint8_t SET_POS_POS_NEG = 64;
static constexpr uint8_t SET_POS_POS_POS = 128;
static constexpr uint8_t CLEAR_NEG_NEG_NEG = 254;
static constexpr uint8_t CLEAR_NEG_NEG_POS = 253;
static constexpr uint8_t CLEAR_NEG_POS_NEG = 251;
static constexpr uint8_t CLEAR_NEG_POS_POS = 247;
static constexpr uint8_t CLEAR_POS_NEG_NEG = 239;
static constexpr uint8_t CLEAR_POS_NEG_POS = 223;
static constexpr uint8_t CLEAR_POS_POS_NEG = 191;
static constexpr uint8_t CLEAR_POS_POS_POS = 127;
static constexpr uint8_t CLEAR_NEG_NEG_ANY =
  (CLEAR_NEG_NEG_NEG & CLEAR_NEG_NEG_POS);
static constexpr uint8_t CLEAR_NEG_POS_ANY = 
  (CLEAR_NEG_POS_NEG & CLEAR_NEG_POS_POS);
static constexpr uint8_t CLEAR_POS_NEG_ANY =
  (CLEAR_POS_NEG_NEG & CLEAR_POS_NEG_POS);
static constexpr uint8_t CLEAR_POS_POS_ANY =
  (CLEAR_POS_POS_NEG & CLEAR_POS_POS_POS);
static constexpr uint8_t CLEAR_NEG_ANY_NEG =
  (CLEAR_NEG_NEG_NEG & CLEAR_NEG_POS_NEG);
static constexpr uint8_t CLEAR_NEG_ANY_POS =
  (CLEAR_NEG_NEG_POS & CLEAR_NEG_POS_POS);
static constexpr uint8_t CLEAR_POS_ANY_NEG =
  (CLEAR_POS_NEG_NEG & CLEAR_POS_POS_NEG);
static constexpr uint8_t CLEAR_POS_ANY_POS =
  (CLEAR_POS_NEG_POS & CLEAR_POS_POS_POS);
static constexpr uint8_t CLEAR_ANY_NEG_NEG =
  (CLEAR_NEG_NEG_NEG & CLEAR_POS_NEG_NEG);
static constexpr uint8_t CLEAR_ANY_NEG_POS =
  (CLEAR_NEG_NEG_POS & CLEAR_POS_NEG_POS);
static constexpr uint8_t CLEAR_ANY_POS_NEG =
  (CLEAR_NEG_POS_NEG & CLEAR_POS_POS_NEG);
static constexpr uint8_t CLEAR_ANY_POS_POS =
  (CLEAR_NEG_POS_POS & CLEAR_POS_POS_POS);
static constexpr uint8_t CLEAR_NEG_ANY_ANY =
  (CLEAR_NEG_NEG_ANY & CLEAR_NEG_POS_ANY);
static constexpr uint8_t CLEAR_POS_ANY_ANY =
  (CLEAR_POS_NEG_ANY & CLEAR_POS_POS_ANY);
static constexpr uint8_t CLEAR_ANY_NEG_ANY =
  (CLEAR_NEG_NEG_ANY & CLEAR_POS_NEG_ANY);
static constexpr uint8_t CLEAR_ANY_POS_ANY =
  (CLEAR_NEG_POS_ANY & CLEAR_POS_POS_ANY);
static constexpr uint8_t CLEAR_ANY_ANY_NEG =
  (CLEAR_NEG_ANY_NEG & CLEAR_POS_ANY_NEG);
static constexpr uint8_t CLEAR_ANY_ANY_POS =
  (CLEAR_NEG_ANY_POS & CLEAR_POS_ANY_POS);

// Aggregate basis state constants for first term (i)
static constexpr uint8_t SET_NEG_ANY_ANY = SET_NEG_NEG_NEG | SET_NEG_NEG_POS | SET_NEG_POS_NEG | SET_NEG_POS_POS;  // 15
static constexpr uint8_t SET_POS_ANY_ANY = SET_POS_NEG_NEG | SET_POS_NEG_POS | SET_POS_POS_NEG | SET_POS_POS_POS;  // 240

// Aggregate basis state constants for second term (j)
static constexpr uint8_t SET_ANY_NEG_ANY = SET_NEG_NEG_NEG | SET_NEG_NEG_POS | SET_POS_NEG_NEG | SET_POS_NEG_POS;  // 51
static constexpr uint8_t SET_ANY_POS_ANY = SET_NEG_POS_NEG | SET_NEG_POS_POS | SET_POS_POS_NEG | SET_POS_POS_POS;  // 204

// Aggregate basis state constants for third term (k)
static constexpr uint8_t SET_ANY_ANY_NEG = SET_NEG_NEG_NEG | SET_NEG_POS_NEG | SET_POS_NEG_NEG | SET_POS_POS_NEG;  // 85
static constexpr uint8_t SET_ANY_ANY_POS = SET_NEG_NEG_POS | SET_NEG_POS_POS | SET_POS_NEG_POS | SET_POS_POS_POS;  // 170

// Additional aggregate basis state constants for specific combinations
// For i,j pair constants
static constexpr uint8_t SET_NEG_NEG_ANY = SET_NEG_NEG_NEG | SET_NEG_NEG_POS;  // 3
static constexpr uint8_t SET_NEG_POS_ANY = SET_NEG_POS_NEG | SET_NEG_POS_POS;  // 12
static constexpr uint8_t SET_POS_NEG_ANY = SET_POS_NEG_NEG | SET_POS_NEG_POS;  // 48
static constexpr uint8_t SET_POS_POS_ANY = SET_POS_POS_NEG | SET_POS_POS_POS;  // 192

// For i,k pair constants
static constexpr uint8_t SET_NEG_ANY_NEG = SET_NEG_NEG_NEG | SET_NEG_POS_NEG;  // 5
static constexpr uint8_t SET_NEG_ANY_POS = SET_NEG_NEG_POS | SET_NEG_POS_POS;  // 10
static constexpr uint8_t SET_POS_ANY_NEG = SET_POS_NEG_NEG | SET_POS_POS_NEG;  // 80
static constexpr uint8_t SET_POS_ANY_POS = SET_POS_NEG_POS | SET_POS_POS_POS;  // 160

// For j,k pair constants
static constexpr uint8_t SET_ANY_NEG_NEG = SET_NEG_NEG_NEG | SET_POS_NEG_NEG;  // 17
static constexpr uint8_t SET_ANY_NEG_POS = SET_NEG_NEG_POS | SET_POS_NEG_POS;  // 34
static constexpr uint8_t SET_ANY_POS_NEG = SET_NEG_POS_NEG | SET_POS_POS_NEG;  // 68
static constexpr uint8_t SET_ANY_POS_POS = SET_NEG_POS_POS | SET_POS_POS_POS;  // 136
static constexpr uint8_t SET_ANY_ANY_ANY = SET_ANY_NEG_NEG | SET_ANY_NEG_POS | SET_ANY_POS_NEG | SET_ANY_POS_POS; // 256

// plug in the positive or negative-ness of the terms and get the
// state value that represents that combination.

static constexpr uint8_t threed_clear_masks[2][2][2] =
  { { { CLEAR_NEG_NEG_NEG, CLEAR_NEG_NEG_POS },
      { CLEAR_NEG_POS_NEG, CLEAR_NEG_POS_POS },},
    { { CLEAR_POS_NEG_NEG, CLEAR_POS_NEG_POS },
      { CLEAR_POS_POS_NEG, CLEAR_POS_POS_POS },},};

static constexpr uint8_t threed_set_masks[2][2][2] =
  { { { SET_NEG_NEG_NEG, SET_NEG_NEG_POS },
      { SET_NEG_POS_NEG, SET_NEG_POS_POS },},
    { { SET_POS_NEG_NEG, SET_POS_NEG_POS },
      { SET_POS_POS_NEG, SET_POS_POS_POS },},};
