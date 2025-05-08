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

#include "pairing.h"
#include <cmath>

/**
 * @brief Maps a pair of indices (i,j) where i < j to a unique single index
 * 
 * This function implements a pairing function that generates a unique index for
 * each unique pair of indices (i,j) where i < j, allowing us to store
 * pair-based data in a flat array. This is particularly useful for
 * storing pair_states in the 3SAT solver algorithm.
 * 
 * The formula used is based on the triangular number formula: (j * (j
 * - 1)) / 2 + i 
 * 
 * For example:
 * - pair2d(0,1) = (1 * (1 - 1)) / 2 + 0 = 0
 * - pair2d(0,2) = (2 * (2 - 1)) / 2 + 0 = 1
 * - pair2d(1,2) = (2 * (2 - 1)) / 2 + 1 = 2
 * - pair2d(0,3) = (3 * (3 - 1)) / 2 + 0 = 3
 * - and so on...
 * 
 * @param i First index (must be less than j)
 * @param j Second index (must be greater than i)
 * @return uint64_t A unique index for the pair
 */
uint64_t pair2d(uint64_t i, uint64_t j) {
    return (j * (j - 1)) / 2 + i;
}

/**
 * @brief Inverse of pair2d - maps a single index back to the original
 * pair (i,j)
 * 
 * This function reverses the pairing function, taking a single index and
 * reconstructing the original pair of indices (i,j) where i < j. This is used
 * when we need to identify which specific pair corresponds to a given index
 * in the flat pair_states array.
 * 
 * The algorithm works as follows:
 * 1. Find j by solving the quadratic equation: (j * (j - 1)) / 2 ≤
 *  index < (j * (j + 1)) / 2 
 * 2. Once j is found, compute i = index - (j * (j - 1)) / 2
 * 
 * @param index The flat array index to convert back to a pair
 * @return std::tuple<uint64_t, uint64_t> The original (i,j) pair
 */
std::tuple<uint64_t, uint64_t> unpair2d(uint64_t index) {
    // Compute discriminant for the quadratic formula
    double discriminant = 1.0 + 8.0 * index;
    
    // Use the quadratic formula to get an approximation of j
    // j ≈ (1 + √(1 + 8*index)) / 2
    double j_approx = (1.0 + std::sqrt(discriminant)) / 2.0;
    uint64_t j = static_cast<uint64_t>(j_approx);
    
    // Refine the estimate of j to ensure it's correct
    // Adjust j downward if needed
    while ((j * (j - 1)) / 2 > index) j--;
    
    // Adjust j upward if needed
    while ((j + 1) * j / 2 <= index) j++;
    
    // Calculate i using the relationship index = (j * (j - 1)) / 2 + i
    uint64_t i = index - (j * (j - 1)) / 2;
    
    return std::make_tuple(i, j);
}

/**
 * @brief Maps a triplet of indices (i,j,k) where i < j < k to a
 * unique single index 
 * 
 * This function is an extension of pair2d for three indices. It
 * generates a unique index for each unique triplet of indices (i,j,k)
 * where i < j < k. This allows us to store triplet-based data (basis
 * states) in a flat array. 
 * 
 * The formula uses a combination of the triangular number and
 * tetrahedral number 
 * formulas: (k * (k - 1) * (k - 2)) / 6 + (j * (j - 1)) / 2 + i
 * 
 * @param i First index (must be less than j)
 * @param j Second index (must be between i and k)
 * @param k Third index (must be greater than j)
 * @return uint64_t A unique index for the triplet
 */
uint64_t pair3d(uint64_t i, uint64_t j, uint64_t k) {
    return (k * (k - 1) * (k - 2)) / 6 + (j * (j - 1)) / 2 + i;
}

/**
 * @brief Inverse of pair3d - maps a single index back to the original
 * triplet (i,j,k) 
 * 
 * This function reverses the triplet pairing function, taking a
 * single index and reconstructing the original triplet of indices
 * (i,j,k) where i < j < k. This is used when we need to identify
 * which specific basis triplet corresponds to a given index in the
 * flat basis_states array. 
 * 
 * The algorithm works as follows:
 * 1. Find k by solving the cubic equation: (k * (k - 1) * (k - 2)) / 6 ≤ index
 * 2. Calculate the remaining index after removing k's contribution
 * 3. Use unpair2d's logic to recover j and i from the remaining index
 * 
 * @param index The flat array index to convert back to a triplet
 * @return std::tuple<uint64_t, uint64_t, uint64_t> The original (i,j,k) triplet
 */
std::tuple<uint64_t, uint64_t, uint64_t> unpair3d(uint64_t index) {
    // Approximate k using the cube root formula
    // This is an approximation for solving (k * (k - 1) * (k - 2)) / 6 = index
    double k_approx = std::cbrt(6.0 * index);
    uint64_t k = static_cast<uint64_t>(k_approx);
    
    // Refine the estimate of k to ensure it's correct
    // Adjust k downward if needed
    while ((k * (k - 1) * (k - 2)) / 6 > index) k--;
    
    // Adjust k upward if needed
    while ((k + 1) * k * (k - 1) / 6 <= index) k++;
    
    // Calculate the remaining index after removing k's contribution
    uint64_t remaining = index - (k * (k - 1) * (k - 2)) / 6;
    
    // Use the same approach as unpair2d to recover j and i
    double discriminant = 1.0 + 8.0 * remaining;
    double j_approx = (1.0 + std::sqrt(discriminant)) / 2.0;
    uint64_t j = static_cast<uint64_t>(j_approx);
    
    // Refine the estimate of j
    while ((j * (j - 1)) / 2 > remaining) j--;
    while ((j + 1) * j / 2 <= remaining) j++;
    
    // Calculate i using the relationship remaining = (j * (j - 1)) / 2 + i
    uint64_t i = remaining - (j * (j - 1)) / 2;
    
    return std::make_tuple(i, j, k);
}

/**
 * @brief Calculates the array size needed to store all pairs (i,j)
 * where i < j < n 
 * 
 * This function computes the number of unique pairs (i,j) where i < j < n.
 * This is equivalent to the binomial coefficient (n choose 2), which equals
 * n! / (2! * (n-2)!) = n * (n - 1) / 2.
 * 
 * This function is used to allocate the proper amount of memory for
 * the pair_states array in the 3SAT solver algorithm.
 * 
 * @param n The upper bound (exclusive) for indices
 * @return uint64_t The number of unique pairs
 */
uint64_t calculate_array_size_2d(uint64_t n) {
    return (n * (n - 1)) / 2;
}

/**
 * @brief Calculates the array size needed to store all triplets
 * (i,j,k) where i < j < k < n 
 * 
 * This function computes the number of unique triplets (i,j,k) where
 * i < j < k < n. 
 * This is equivalent to the binomial coefficient (n choose 3), which equals
 * n! / (3! * (n-3)!) = n * (n - 1) * (n - 2) / 6.
 * 
 * This function is used to allocate the proper amount of memory for
 * the basis_states array in the 3SAT solver algorithm.
 * 
 * @param n The upper bound (exclusive) for indices
 * @return uint64_t The number of unique triplets
 */
uint64_t calculate_array_size_3d(uint64_t n) {
    return (n * (n - 1) * (n - 2)) / 6;
}
