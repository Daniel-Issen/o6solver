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

#pragma once
#include <cstdint>
#include <vector>
#include <tuple>
#include "constants.h"
#include "pairing.h"

// Result structure for state updates
struct UpdateResult {
    bool changed;     // True if any state was changed
    bool has_zero;    // True if any state equals zero
    
    // Constructor
    UpdateResult(bool c = false, bool z = false) : changed(c), has_zero(z) {}
    
    // Convenience function to check if either condition is true
    bool any() const { return changed || has_zero; }
};

std::string term_state_str(uint8_t state);
std::string pair_state_str(uint8_t state);
std::string basis_state_str(uint8_t state);

// update pair and term states to maintain consistency
UpdateResult update_pair_states(Index i, Index j,
				std::vector<uint8_t>& term_states,
				std::vector<uint8_t>& pair_states);

// Update basis, pair, and term states to maintain consistency
// Returns detailed result about what happened during update
UpdateResult update_basis_states(Index i, Index j, Index k,
				 Index basis_idx,
				 std::vector<uint8_t>& term_states,
				 std::vector<uint8_t>& pair_states,
				 std::vector<uint8_t>& basis_states);

UpdateResult ensure_basis_consistency(
    Index basis1_idx, Index basis2_idx, 
    std::vector<uint8_t>& term_states,
    std::vector<uint8_t>& pair_states,
    std::vector<uint8_t>& basis_states);

// Ensure consistency across all bases in the system
bool ensure_global_consistency(std::vector<uint8_t>& term_states,
			       std::vector<uint8_t>& pair_states,
			       std::vector<uint8_t>& basis_states,
			       bool& has_contradiction,
			       Index starting_basis_pair,
			       Index ending_basis_pair);

bool parallel_ensure_global_consistency
(std::vector<uint8_t>& term_states,
 std::vector<uint8_t>& pair_states,
 std::vector<uint8_t>& basis_states,
 bool& has_contradiction,
 Index starting_basis_pair,
 Index ending_basis_pair,
 int num_workers);


