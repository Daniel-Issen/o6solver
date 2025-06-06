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

#include <vector>
#include <string>
#include <cstdint>
#include "file_parser.h"

// Directly apply CNF constraints without creating unnecessary dummy variables
bool apply_constraints(const std::vector<std::vector<Literal>>& cnf_clauses,
		       int& num_vars,
		       std::vector<uint8_t>& term_states,
		       std::vector<uint8_t>& pair_states,
		       std::vector<uint8_t>& basis_states);

// Check satisfiability using the optimized approach
bool check_satisfiability
(int num_workers,
 const std::vector<std::vector<Literal>>& cnf_clauses, 
 int num_vars, 
 bool find_solution = false,
 const std::string& solution_file = "");

// Cross-level consistency checking
bool ensure_cross_level_consistency(std::vector<uint8_t>& term_states,
                                   std::vector<uint8_t>& pair_states,
                                   std::vector<uint8_t>& basis_states);

