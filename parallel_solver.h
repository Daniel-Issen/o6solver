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
#include <thread>
#include <future>
#include <cstdint>
#include "basis_consistency.h"
#include "pairing.h"
#include "file_parser.h"

// Structure to hold work segment boundaries
struct WorkSegment {
  uint64_t starting_basis_pair;
  uint64_t ending_basis_pair;
};

// Structure to hold worker results
struct WorkerResult {
  uint64_t updates;
  bool has_contradiction;
  std::vector<uint8_t> term_states;
  std::vector<uint8_t> pair_states;
  std::vector<uint8_t> basis_states;
};

// Function to divide work among workers
std::vector<WorkSegment> divide_work(uint64_t n, int num_workers);

// Worker function that processes a segment
WorkerResult process_segment(const WorkSegment& segment,
			     const std::vector<uint8_t>& term_states,
			     const std::vector<uint8_t>& pair_states,
			     const std::vector<uint8_t>& basis_states);

// Main parallel solver function
bool parallel_check_satisfiability
(const std::vector<std::vector<Literal>>& cnf_clauses,
 int num_vars,
 int num_workers,
 bool find_solution = false,
 const std::string& solution_file = "");

