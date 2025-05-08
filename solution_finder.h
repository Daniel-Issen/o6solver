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
#include <set>
#include <string>
#include "basis_consistency.h"
#include "file_parser.h"

// Structure to hold a 3SAT solution
struct SATSolution {
    std::vector<int8_t> assignments; // -1 = NEG, 1 = POS, 0 = Unassigned
};

// Function to validate that a solution satisfies all clauses in the formula
bool validate_solution(const SATSolution& solution, 
                       const std::vector<std::vector<Literal>>& cnf_clauses);

// Function to determine a solution from the current states
SATSolution determine_solution(
    std::vector<uint8_t>& basis_states,
    std::vector<uint8_t>& pair_states,
    std::vector<uint8_t>& term_states,
    uint64_t n);

// Helper function to save solution to a file
bool save_solution_to_file(const SATSolution& solution,
			   const std::string& filename);

// Helper function to print solution to console
void print_solution(const SATSolution& solution);
